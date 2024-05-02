import pytest
from cobra import Model
from eflux.eflux2 import add_slack_variables_to_model


@pytest.fixture
def model():
    """Create cobra model."""
    return Model()


# TODO: Check these tests for relavance to eflux (Note: add_slack_variables_to_model is still under development)
def test_add_slack_variables_to_model(model, upper_bounds):
    """Test add_slack_variables_to_model function."""
    # Test case 1: model and upper_bounds are valid inputs
    relaxed_model = add_slack_variables_to_model(model, upper_bounds)
    assert len(relaxed_model.variables) == len(upper_bounds)

    # Test case 2: model is None
    with pytest.raises(TypeError):
        add_slack_variables_to_model(None, upper_bounds)

    # Test case 3: upper_bounds is None
    with pytest.raises(TypeError):
        add_slack_variables_to_model(model, None)

    # Test case 4: model and upper_bounds are empty
    relaxed_model = add_slack_variables_to_model(model, {})
    assert len(relaxed_model.variables) == 0
