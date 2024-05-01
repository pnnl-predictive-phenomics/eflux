"""Tests for eflux.utils."""
import cobra
import pandas as pd
import pytest
from eflux.utils import get_flux_bounds, get_gpr_dict


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

@pytest.fixture(
    name="cobra_model_with_one_gene_reaction_rule",
)
def add_genes_to_r2(cobra_model):
    r2 = cobra_model.reactions.get_by_id("r2")
    r2.gene_reaction_rule = "(gene1 and gene2) or (gene3 and gene4)"
    return cobra_model

@pytest.fixture(
    name="cobra_model_with_multiple_gene_reaction_rule",
)
def add_genes_to_r3(cobra_model_with_one_gene_reaction_rule):
    r3 = cobra_model_with_one_gene_reaction_rule.reactions.get_by_id("r3")
    r3.gene_reaction_rule = "(gene5 and gene6) or (gene7 and gene8)"
    return cobra_model_with_one_gene_reaction_rule

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


def test_empty_model():
    model = cobra.Model("test_model")
    gpr_dict = get_gpr_dict(model)
    assert gpr_dict == {}


def test_model_with_no_gene_reaction_rule(cobra_model):
    gpr_dict = get_gpr_dict(cobra_model)
    assert gpr_dict == {}


def test_model_with_gene_reaction_rule(cobra_model_with_one_gene_reaction_rule):
    gpr_dict = get_gpr_dict(cobra_model_with_one_gene_reaction_rule)
    assert gpr_dict == {'r2': {frozenset(["gene1", "gene2"]), frozenset(["gene3", "gene4"])}}


def test_model_with_multiple_gene_reaction_rules(cobra_model_with_multiple_gene_reaction_rule):
    gpr_dict = get_gpr_dict(cobra_model_with_multiple_gene_reaction_rule)
    assert gpr_dict == {'r2': {frozenset(["gene1", "gene2"]), frozenset(["gene3", "gene4"])},
                        'r3': {frozenset(["gene5", "gene6"]), frozenset(["gene7", "gene8"])}}


