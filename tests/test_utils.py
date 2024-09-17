"""Tests for eflux.utils."""

import numpy as np
import pandas as pd
import pytest
from cobra.core.model import Model
from eflux.utils import (
    convert_transcriptomics_to_enzyme_activity,
    gene_expression_to_enzyme_activity,
    get_gpr_dict,
    get_max_flux_bounds,
    load_model_from_path,
)


def test_get_max_flux_bounds(cobra_model):
    """Test flux bounds for r1 and r2."""
    flux_bounds = get_max_flux_bounds(cobra_model, ["r1", "r4"])

    assert flux_bounds["r2"] == 5
    assert flux_bounds["r3"] == 5


def test_get_max_flux_bounds_with_zero_threshold(cobra_model):
    """Test flux bounds for r1 and r2 with a zero threshold defined."""
    # Get flux bounds for r1 and r2 with a zero threshold
    flux_bounds = get_max_flux_bounds(cobra_model, ["r1", "r4"], precision=1)

    # Check that flux bounds are set correctly
    assert flux_bounds["r2"] == 5
    assert flux_bounds["r3"] == 5

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
    assert gpr_dict == {
        r2: {frozenset(["gene1", "gene2"])},
        r3: {frozenset(["gene5", "gene6"]), frozenset(["gene7", "gene8"])},
    }


def test_gene_expression_to_enzyme_activity(cobra_model_2, expression):
    """Test gene_expression_to_enzyme_activity function in utils.py."""
    # Test enzyme_activity for reaction with no genes.
    r1 = cobra_model_2.reactions.get_by_id("r1")
    gpr0 = get_gpr_dict(cobra_model_2)
    result = gene_expression_to_enzyme_activity(cobra_model_2, gpr0, expression)
    assert np.isnan(result[r1])

    # Test enzyme_activity for reaction with one gene.
    r1.gene_reaction_rule = "gene3"
    gpr1 = get_gpr_dict(cobra_model_2)
    result = gene_expression_to_enzyme_activity(cobra_model_2, gpr1, expression)
    assert result[r1] == expression["gene3"]

    # Test enzyme_activity for reaction with multiple genes.
    r2 = cobra_model_2.reactions.get_by_id("r2")
    gpr2 = get_gpr_dict(cobra_model_2)
    expression["gene1"] = 0.0
    expression["gene2"] = 0.0
    result = gene_expression_to_enzyme_activity(cobra_model_2, gpr2, expression)
    assert result[r2] == 0.0

    # Test enzyme_activity for reaction with isozyme.
    r3 = cobra_model_2.reactions.get_by_id("r3")
    gpr3 = get_gpr_dict(cobra_model_2)
    expression["gene7"] = 0.0
    expression["gene8"] = 0.0
    result = gene_expression_to_enzyme_activity(cobra_model_2, gpr3, expression)
    assert result[r3] == np.min([expression["gene5"], expression["gene6"]])

    # Test enzyme_activity for reaction with unobserved gene.
    r4 = cobra_model_2.reactions.get_by_id("r4")
    r4.gene_reaction_rule = "gene9"
    gpr4 = get_gpr_dict(cobra_model_2)
    result = gene_expression_to_enzyme_activity(cobra_model_2, gpr4, expression)
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


def test_convert_transcriptomics_to_enzyme_activity(
    cobra_model_2, input_transcriptomics, expected_enzyme_activity
):
    """Test convert_transcriptomics_to_enzyme_activity function in utils.py."""
    # Test convert_transcriptomics_to_enzyme_activity for empty data and empty model
    result = convert_transcriptomics_to_enzyme_activity(pd.DataFrame(), Model())
    assert result.empty

    # Test convert_transcriptomics_to_enzyme_activity for non-empty data and empty model
    result = convert_transcriptomics_to_enzyme_activity(input_transcriptomics, Model())
    assert result.empty

    # Test convert_transcriptomics_to_enzyme_activity for empty data and non-empty model
    result = convert_transcriptomics_to_enzyme_activity(pd.DataFrame(), cobra_model_2)
    assert result.empty

    # Test convert_transcriptomics_to_enzyme_activity for non-empty data and non-empty model
    result = convert_transcriptomics_to_enzyme_activity(input_transcriptomics, cobra_model_2)
    assert result.shape == (4, 2)
    assert result.equals(expected_enzyme_activity)


def test_load_model_success(valid_model_paths):
    """Test successful loading of models from different formats."""
    for _, path in valid_model_paths.items():
        model = load_model_from_path(path)
        assert isinstance(model, Model)


def test_load_model_invalid_extension():
    """Test handling of unsupported file extensions."""
    with pytest.raises(ValueError) as exc_info:
        load_model_from_path("invalid_model.txt")
    assert "Unsupported model file extension" in str(exc_info.value)
