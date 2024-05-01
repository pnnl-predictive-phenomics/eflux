"""Tests for eflux.utils."""
import cobra
import pandas as pd
import pytest
from eflux.utils import get_flux_bounds


@pytest.fixture(
    name="cobra_model",
)
def get_cobra_model():
    """Fixture for testing a toy cobra model."""
    model = cobra.Model("test_model")
    r1 = cobra.Reaction("r1")
    r2 = cobra.Reaction("r2")
    r3 = cobra.Reaction("r3")
    r4 = cobra.Reaction("r4")

    model.add_metabolites([cobra.Metabolite("m1"), cobra.Metabolite("m2"), cobra.Metabolite("m3")])
    r1.add_metabolites({model.metabolites.m1: 1})
    r2.add_metabolites({model.metabolites.m1: -1, model.metabolites.m2: 1})
    r3.add_metabolites({model.metabolites.m2: -1, model.metabolites.m3: 1})
    r4.add_metabolites({model.metabolites.m3: -1})

    # Set flux bounds for r1 and r2
    r2.bounds = (0.01, 10)
    r3.bounds = (0, 5)

    model.add_reactions([r1, r2, r3, r4])
    return model


def test_get_flux_bounds(cobra_model):
    # Get flux bounds for r1 and r2
    model, flux_bounds = get_flux_bounds(cobra_model, ["r1", "r4"])

    # Check that flux bounds are set correctly
    assert flux_bounds.loc["r2", "minimum"] == 0.01
    assert flux_bounds.loc["r2", "maximum"] == 5
    assert flux_bounds.loc["r3", "minimum"] == 0.01
    assert flux_bounds.loc["r3", "maximum"] == 5

def test_get_flux_bounds_with_zero_threshold(cobra_model):
    # Get flux bounds for r1 and r2 with a zero threshold
    model, flux_bounds = get_flux_bounds(cobra_model, ["r1", "r4"], zero_threshold=0.1)

    # Check that flux bounds are set correctly
    assert flux_bounds.loc["r2", "minimum"] == 0
    assert flux_bounds.loc["r2", "maximum"] == 5
    assert flux_bounds.loc["r3", "minimum"] == 0
    assert flux_bounds.loc["r3", "maximum"] == 5

    # Check that flux bounds are set correctly when zero_threshold is set to 0