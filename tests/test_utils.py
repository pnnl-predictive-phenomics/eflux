"""Tests for eflux.utils."""
import numpy as np
import pandas as pd
import pytest
from cobra.core import Gene, Metabolite, Reaction
from cobra.core.model import Model
from eflux.utils import (
    convert_transcriptomics_to_enzyme_activity,
    gene_expression_to_enzyme_activity,
    get_flux_bounds,
    get_gpr_dict,
)
from fixtures import (
    add_genes_to_r2,
    add_genes_to_r3,
    expression,
    get_cobra_model,
    transcriptomics_data,
)


def test_get_flux_bounds(cobra_model):
    """Test flux bounds for r1 and r2."""
    model, flux_bounds = get_flux_bounds(cobra_model, ["r1", "r4"])

    assert flux_bounds.loc["r2", "minimum"] == 0.01
    assert flux_bounds.loc["r2", "maximum"] == 5
    assert flux_bounds.loc["r3", "minimum"] == 0.01
    assert flux_bounds.loc["r3", "maximum"] == 5


def test_get_flux_bounds_with_zero_threshold(cobra_model):
    """Test flux bounds for r1 and r2 with a zero threshold defined."""
    # Get flux bounds for r1 and r2 with a zero threshold
    model, flux_bounds = get_flux_bounds(cobra_model, ["r1", "r4"], zero_threshold=0.1)

    # Check that flux bounds are set correctly
    assert flux_bounds.loc["r2", "minimum"] == 0
    assert flux_bounds.loc["r2", "maximum"] == 5
    assert flux_bounds.loc["r3", "minimum"] == 0
    assert flux_bounds.loc["r3", "maximum"] == 5

    # Check that flux bounds are set correctly when zero_threshold is set to 0


def test_gpr_dict_for_empty_model():
    """Test get_gpr_dict for an empty model."""
    model = Model("test_model")
    gpr_dict = get_gpr_dict(model)
    assert gpr_dict == {}


def test_gpr_dict_for_model_with_no_gene_reaction_rule(cobra_model):
    """Test get_gpr_dict for a model with no gene reaction rule."""
    gpr_dict = get_gpr_dict(cobra_model)
    assert gpr_dict == {}


def test_gpr_dict_for_model_with_gene_reaction_rule(cobra_model_1):
    """Test get_gpr_dict for a model with a single gene reaction rule."""
    gpr_dict = get_gpr_dict(cobra_model_1)
    r2 = cobra_model_1.reactions.get_by_id("r2")
    assert gpr_dict == {r2: {frozenset(["gene1", "gene2"])}}


def test_gpr_dict_for_model_with_multiple_gene_reaction_rules(cobra_model_2):
    """Test get_gpr_dict for a model with multiple gene reaction rules (i.e. isozyme)."""
    gpr_dict = get_gpr_dict(cobra_model_2)
    r2 = cobra_model_2.reactions.get_by_id("r2")
    r3 = cobra_model_2.reactions.get_by_id("r3")
    assert gpr_dict == {r2: {frozenset(["gene1", "gene2"])},
                        r3: {frozenset(["gene5", "gene6"]), frozenset(["gene7", "gene8"])}}


def test_enzyme_activity_for_no_genes(cobra_model_2, expression):
    """Test enzyme_activity for reaction with no genes."""
    r1 = cobra_model_2.reactions.get_by_id("r1")
    gpr = get_gpr_dict(cobra_model_2)
    result = gene_expression_to_enzyme_activity(cobra_model_2, gpr, expression)
    assert np.isnan(result[r1])

def test_enzyme_activity_for_one_gene(cobra_model_2, expression):
    """Test enzyme_activity for reaction with one gene."""
    r1 = cobra_model_2.reactions.get_by_id("r1")
    r1.gene_reaction_rule = "gene3"
    gpr = get_gpr_dict(cobra_model_2)
    result = gene_expression_to_enzyme_activity(cobra_model_2, gpr, expression)
    assert result[r1] == expression["gene3"]

def test_enzyme_activity_for_multiple_genes(cobra_model_2, expression):
    """Test enzyme_activity for reaction with multiple genes."""
    r2 = cobra_model_2.reactions.get_by_id("r2")
    gpr = get_gpr_dict(cobra_model_2)
    expression["gene1"] = 0.0
    expression["gene2"] = 0.0
    result = gene_expression_to_enzyme_activity(cobra_model_2, gpr, expression)
    assert result[r2] == 0.0

def test_enzyme_activity_multiple_isozymes(cobra_model_2, expression):
    """Test enzyme_activity for reaction with isozyme."""
    r3 = cobra_model_2.reactions.get_by_id("r3")
    gpr = get_gpr_dict(cobra_model_2)
    expression["gene7"] = 0.0
    expression["gene8"] = 0.0
    result = gene_expression_to_enzyme_activity(cobra_model_2, gpr, expression)
    assert result[r3] == np.min([expression["gene5"], expression["gene6"]])

def test_enzyme_activity_unobserved_gene(cobra_model_2, expression):
    """Test enzyme_activity for reaction with unobserved gene."""
    r4 = cobra_model_2.reactions.get_by_id("r4")
    r4.gene_reaction_rule = "gene9"
    gpr = get_gpr_dict(cobra_model_2)
    result = gene_expression_to_enzyme_activity(cobra_model_2, gpr, expression)
    assert result[r4] == np.inf




# @pytest.fixture
# def model():
#     model = Model()
#     model.add_reaction(cobra.Reaction('rxn1'))
#     model.add_reaction(cobra.Reaction('rxn2'))
#     model.reactions['rxn1'].gene_reaction_rule = 'gene1 OR gene2'
#     model.reactions['rxn2'].gene_reaction_rule = 'gene3 AND gene4'
#     model.add_genes([cobra.Gene('gene1'), cobra.Gene('gene2'), cobra.Gene('gene3'), cobra.Gene('gene4')])
#     return model

def test_empty_data_empty_model(transcriptomics, cobra_model):
    """Test convert_transcriptomics_to_enzyme_activity for empty data and empty model."""
    result = convert_transcriptomics_to_enzyme_activity(pd.DataFrame(), Model())
    assert result.empty


def test_non_empty_data_empty_model(transcriptomics, cobra_model):
    """Test convert_transcriptomics_to_enzyme_activity for non-empty data and empty model."""
    result = convert_transcriptomics_to_enzyme_activity(transcriptomics, Model())
    assert result.empty


def test_empty_data_non_empty_model(transcriptomics, cobra_model_2):
    """Test convert_transcriptomics_to_enzyme_activity for empty data and non-empty model."""
    result = convert_transcriptomics_to_enzyme_activity(pd.DataFrame(), cobra_model_2)
    assert result.empty


def test_non_empty_data_non_empty_model(transcriptomics, cobra_model_2):
    """Test convert_transcriptomics_to_enzyme_activity for non-empty data and non-empty model."""
    result = convert_transcriptomics_to_enzyme_activity(transcriptomics, cobra_model_2)
    expected_df = pd.DataFrame({'Reaction_ID': ['r1', 'r2', 'r3', 'r4'],
                                'strain1': [np.nan, 1.0, np.inf, np.nan],
                                'strain2': [np.nan, 4.0, np.inf, np.nan],
                                })
    expected_df = expected_df.set_index('Reaction_ID')
    assert result.shape == (4, 2)
    assert result.equals(expected_df)

