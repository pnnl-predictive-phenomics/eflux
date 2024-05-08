"""Tests for eflux functions."""

import pytest
from cobra import exceptions
from eflux.eflux2 import add_slack_variables_to_model, get_normalized_condition, get_upper_bounds


def test_add_slack_variables_to_model(min_uptake_model, infeasible_upper_bounds, expected_fluxes):
    """Test add_slack_variables_to_model function."""
    # Test case 1: model is None
    with pytest.raises(TypeError):
        add_slack_variables_to_model(None, infeasible_upper_bounds)

    # Test case 2: upper_bounds is None
    with pytest.raises(TypeError):
        add_slack_variables_to_model(min_uptake_model, None)

    # Test case 3: model and upper_bounds are empty
    relaxed_model = add_slack_variables_to_model(min_uptake_model, {})
    assert len([v for v in relaxed_model.variables if "SLACK" in v.name]) == 0

    # Test case 4: model and upper_bounds are valid inputs
    relaxed_model = add_slack_variables_to_model(min_uptake_model, infeasible_upper_bounds)
    actual_fluxes = relaxed_model.optimize().fluxes.to_dict()
    slack_vars = [v for v in relaxed_model.variables if "SLACK" in v.name]
    assert len(slack_vars) == len(infeasible_upper_bounds)
    assert slack_vars[0].primal == 1.0
    assert actual_fluxes == expected_fluxes


def test_infeasible_model(infeasible_model):
    """Test infeasible_model function."""
    with pytest.raises(exceptions.OptimizationError):
        infeasible_model.optimize(raise_error=True)


def test_get_normalized_condition(input_enzyme_activity, good_ref_col, good_target_col):
    """Test get_normalized_condition function."""
    # Test get_normalized_condition with bad reference and target column names.
    with pytest.raises(KeyError):
        get_normalized_condition(df=input_enzyme_activity, ref_col='bad_ref_col', target_col='bad_target_col')
    with pytest.raises(KeyError):
        get_normalized_condition(input_enzyme_activity, ref_col=good_ref_col, target_col='bad_target_col')
    with pytest.raises(KeyError):
        get_normalized_condition(input_enzyme_activity, ref_col='bad_ref_col', target_col=good_target_col)

    # Test get_normalized_condition with zero entries in ref_col
    with pytest.raises(ValueError):
        get_normalized_condition(input_enzyme_activity, ref_col='ref_col_with_zero', target_col=good_target_col)

    # Test get_normalized_condition with inf or nan entries in ref_col."""
    with pytest.raises(ValueError):
        get_normalized_condition(input_enzyme_activity, ref_col='ref_col_with_inf', target_col=good_target_col)
    with pytest.raises(ValueError):
        get_normalized_condition(input_enzyme_activity, ref_col='ref_col_with_nan', target_col=good_target_col)

    # Test get_normalized_condition with inf or nan entries in target_col.
    with pytest.raises(ValueError):
        get_normalized_condition(input_enzyme_activity, ref_col=good_ref_col, target_col='target_col_with_inf')
    with pytest.raises(ValueError):
        get_normalized_condition(input_enzyme_activity, ref_col=good_ref_col, target_col='target_col_with_nan')

    # Test get_normalized_condition with normal case.
    result = get_normalized_condition(input_enzyme_activity, ref_col=good_ref_col, target_col=good_target_col)
    assert result == {'r1': 2.0, 'r2': 2.0, 'r3': 2.0, 'r4': 2.0}


def test_get_upper_bounds(cobra_model, input_normalized_enzyme_activity, expected_dict_from_get_enzyme_bounds):
    """Test get_upper_bounds function."""
    # Test get_upper_bounds with empty model input.
    with pytest.raises(AttributeError):
        get_upper_bounds(model=None, scaling_factors=input_normalized_enzyme_activity)

    # Test get_upper_bounds with empty norm_enzyme_activity input.
    enzyme_bounds = get_upper_bounds(model=cobra_model, scaling_factors={})
    assert enzyme_bounds == {}

    # Test get_upper_bounds with valid inputs.
    enzyme_bounds = get_upper_bounds(model=cobra_model, scaling_factors=input_normalized_enzyme_activity)
    assert isinstance(enzyme_bounds, dict)  # is this needed?
    assert enzyme_bounds == expected_dict_from_get_enzyme_bounds
