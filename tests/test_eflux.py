import pytest
from cobra import Model, exceptions
from eflux.eflux2 import add_slack_variables_to_model


# TODO: Check these tests for relavance to eflux (Note: add_slack_variables_to_model is still under development)
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
