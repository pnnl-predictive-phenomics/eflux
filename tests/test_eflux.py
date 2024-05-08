import pytest
from cobra import Model, exceptions
from eflux.eflux2 import add_slack_variables_to_model, get_enzyme_bounds


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


def test_empty_model_input(input_enzyme_activity):
    """Test get_enzyme_bounds with empty model input."""
    enzyme_bounds = get_enzyme_bounds(model={}, norm_enzyme_activity=input_enzyme_activity)
    assert enzyme_bounds == {}


def test_empty_norm_enzyme_activity_input(cobra_model):
    """Test get_enzyme_bounds with empty norm_enzyme_activity input."""
    enzyme_bounds = get_enzyme_bounds(model=cobra_model, norm_enzyme_activity={})
    assert enzyme_bounds == {}


def test_valid_inputs(cobra_model, input_enzyme_activity, expected_dict_from_get_enzyme_bounds):
    """Test get_enzyme_bounds with valid inputs."""
    enzyme_bounds = get_enzyme_bounds(model=cobra_model, norm_enzyme_activity=input_enzyme_activity)
    assert isinstance(enzyme_bounds, dict)  # is this needed?
    assert enzyme_bounds == expected_dict_from_get_enzyme_bounds
