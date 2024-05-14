"""Tests for eflux functions."""

import pytest
from cobra import exceptions
from eflux.eflux2 import (
    add_slack_variables_to_model,
    get_condition_specific_upper_bounds,
    get_normalized_condition,  # run_condition_specific_eflux,
)


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
        get_normalized_condition(
            df=input_enzyme_activity, ref_col="bad_ref_col", target_col="bad_target_col"
        )
    with pytest.raises(KeyError):
        get_normalized_condition(
            input_enzyme_activity, ref_col=good_ref_col, target_col="bad_target_col"
        )
    with pytest.raises(KeyError):
        get_normalized_condition(
            input_enzyme_activity, ref_col="bad_ref_col", target_col=good_target_col
        )

    # Test get_normalized_condition with zero entries in ref_col
    with pytest.raises(ValueError):
        get_normalized_condition(
            input_enzyme_activity, ref_col="ref_col_with_zero", target_col=good_target_col
        )

    # Test get_normalized_condition with inf or nan entries in ref_col."""
    with pytest.raises(ValueError):
        get_normalized_condition(
            input_enzyme_activity, ref_col="ref_col_with_inf", target_col=good_target_col
        )
    with pytest.raises(ValueError):
        get_normalized_condition(
            input_enzyme_activity, ref_col="ref_col_with_nan", target_col=good_target_col
        )

    # Test get_normalized_condition with inf or nan entries in target_col.
    with pytest.raises(ValueError):
        get_normalized_condition(
            input_enzyme_activity, ref_col=good_ref_col, target_col="target_col_with_inf"
        )
    with pytest.raises(ValueError):
        get_normalized_condition(
            input_enzyme_activity, ref_col=good_ref_col, target_col="target_col_with_nan"
        )

    # Test get_normalized_condition with normal case.
    result = get_normalized_condition(
        input_enzyme_activity, ref_col=good_ref_col, target_col=good_target_col
    )
    assert result == {"r1": 2.0, "r2": 2.0, "r3": 2.0, "r4": 2.0}


def test_get_condition_specific_upper_bounds(
    input_upper_bounds, input_normalized_enzyme_activity, expected_dict_from_get_enzyme_bounds
):
    """Test get_condition_specific_upper_bounds function."""
    # Test get_upper_bounds with empty fva_upper_bounds input.
    with pytest.raises(AttributeError):
        get_condition_specific_upper_bounds(
            fva_upper_bounds=None, scaling_factors=input_normalized_enzyme_activity
        )

    # Test get_upper_bounds with empty norm_enzyme_activity input.
    enzyme_bounds = get_condition_specific_upper_bounds(
        fva_upper_bounds=input_upper_bounds, scaling_factors={}
    )
    assert enzyme_bounds == {}

    # Test get_upper_bounds with valid inputs.
    enzyme_bounds = get_condition_specific_upper_bounds(
        fva_upper_bounds=input_upper_bounds, scaling_factors=input_normalized_enzyme_activity
    )
    assert isinstance(enzyme_bounds, dict)  # is this needed?
    assert enzyme_bounds == expected_dict_from_get_enzyme_bounds


# @pytest.fixture
# def external_fluxes():
#     return pd.DataFrame({'rxn1': [1.0], 'rxn2': [2.0]})


# @pytest.fixture
# def enzyme_activity():
#     return pd.DataFrame({'cond1': [1.0], 'cond2': [2.0]})


# def test_run_condition_specific_eflux(cobra_model_2, external_fluxes, enzyme_activity):
#     """Test run_condition_specific_eflux function."""
#     # Test case 1: model is None
#     with pytest.raises(TypeError):
#         run_condition_specific_eflux(None, 'growth', 'product', external_fluxes, enzyme_activity, 'cond1', 'cond2')

#     # Test case 2: growth_rxn_id is None
#     with pytest.raises(TypeError):
#         run_condition_specific_eflux(cobra_model_2, None, 'product', external_fluxes, enzyme_activity, 'cond1', 'cond2')

#     # Test case 3: product_rxn_id is None
#     with pytest.raises(TypeError):
#         run_condition_specific_eflux(cobra_model_2, 'growth', None, external_fluxes, enzyme_activity, 'cond1', 'cond2')

#     # Test case 4: external_fluxes is None
#     with pytest.raises(TypeError):
#         run_condition_specific_eflux(cobra_model_2, 'growth', 'product', None, enzyme_activity, 'cond1', 'cond2')

#     # Test case 5: enzyme_activity is None
#     with pytest.raises(TypeError):
#         run_condition_specific_eflux(cobra_model_2, 'growth', 'product', external_fluxes, None, 'cond1', 'cond2')

#     # Test case 6: ref_cond is None
#     with pytest.raises(TypeError):
#         run_condition_specific_eflux(cobra_model_2, 'growth', 'product', external_fluxes, enzyme_activity, None, 'cond2')

#     # Test case 7: target_cond is None
#     with pytest.raises(TypeError):
#         run_condition_specific_eflux(cobra_model_2, 'growth', 'product', external_fluxes, enzyme_activity, 'cond1', None)

#     # Test case 8: All inputs are valid
#     assert run_condition_specific_eflux(cobra_model_2, 'growth', 'product', external_fluxes, enzyme_activity, 'cond1', 'cond2') == {}


# Inputs to test eflux
# 1. cobra model

# run_eflux_condition_specific(
#         model: cobra.Model, growth_rxn_id: str, product_rxn_id: str,
#         external_fluxes: pd.DataFrame, enzyme_activity: pd.DataFrame,
#         ref_cond: str, target_cond: str) -> dict[str, float]

# Test data
# Flow for S. elongatus FVA
# 1- Optimize for growth-rate
# 2- Constrain growth-rate at 80% of optimized level --> Optimize for surose production ==> max level to be used as reference flux
