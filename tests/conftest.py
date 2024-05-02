"""Fixtures for eflux tests."""

import pandas as pd
import pytest
from cobra.core import Metabolite, Reaction
from cobra.core.model import Model


@pytest.fixture(
    name="cobra_model",
)
def get_cobra_model():
    """Fixture for testing a toy cobra model."""
    model = Model("test_model")
    r1 = Reaction("r1")
    r2 = Reaction("r2")
    r3 = Reaction("r3")
    r4 = Reaction("r4")

    model.add_metabolites([Metabolite("m1"), Metabolite("m2"), Metabolite("m3")])
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
    name="cobra_model_1",
)
def add_genes_to_r2(cobra_model):
    """Add gene1 and gene2 to Reaction r2."""
    r2 = cobra_model.reactions.get_by_id("r2")
    r2.gene_reaction_rule = "gene1 and gene2"
    return cobra_model


@pytest.fixture(
    name="cobra_model_2",
)
def add_genes_to_r3(cobra_model_1):
    """Add (gene5 and gene6) or (gene7 and gene8) to Reaction r3."""
    r3 = cobra_model_1.reactions.get_by_id("r3")
    r3.gene_reaction_rule = "(gene5 and gene6) or (gene7 and gene8)"
    return cobra_model_1


@pytest.fixture(name="expression")
def expression():
    """Fixture for testing gene expression."""
    return {
        "gene1": 1.0,
        "gene2": 2.0,
        "gene3": 3.0,
        "gene4": 4.0,
        "gene5": 5.0,
        "gene6": 6.0,
        "gene7": 7.0,
        "gene8": 8.0,
    }


@pytest.fixture(name="transcriptomics")
def transcriptomics_data():
    """Fixture for testing transcriptomics data."""
    return pd.DataFrame(
        {"strain1": [1, 2, 3, 4, 5], "strain2": [5, 4, 3, 2, 1]},
        index=["gene1", "gene2", "gene3", "gene5", "gene6"],
    )


@pytest.fixture
def upper_bounds():
    """Fixture to set upper bound on reaction."""
    return {"reaction1": 10.0, "reaction2": 5.0}
