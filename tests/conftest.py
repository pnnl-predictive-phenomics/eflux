"""Fixtures for eflux tests."""

import numpy as np
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


@pytest.fixture(
        name="expected_fluxes",
)
def expected_fluxes():
    """Fixture to set expected fluxes."""
    return {'r1': 4.0, 'r2': 4.0, 'r3': 4.0, 'r4': 4.0}


@pytest.fixture(
        name="model_with_objective",
)
def model_with_objective(cobra_model):
    """Fixture to set reaction r4 as objective."""
    model = cobra_model.copy()
    model.objective = model.reactions.get_by_id('r4')
    return model


@pytest.fixture(
        name="min_uptake_model",
)
def min_uptake_model(model_with_objective):
    """Fixture to set postive lower bound on reaction r1."""
    model = model_with_objective.copy()
    model.reactions.r1.lower_bound = 4.0
    return model


@pytest.fixture(
        name="infeasible_upper_bounds",
)
def infeasible_upper_bounds():
    """Fixture to set infeasible_upper_bounds."""
    return {'r3': 3.0}


@pytest.fixture(
        name="infeasible_model",
)
def infeasible_model(min_uptake_model, infeasible_upper_bounds):
    """Fixture to set infeasible model."""
    model = min_uptake_model.copy()
    for r, b in infeasible_upper_bounds.items():
        model.reactions.get_by_id(r).upper_bound = b
    return model


@pytest.fixture(
        name="input_enzyme_activity",
)
def input_enzyme_activity():
    """Fixture enzyme activity input."""
    return pd.DataFrame({'reference_cond': [1.0, 2.0, 3.0, 4.0],
                         'target_cond': [2.0, 4.0, 6.0, 8.0],
                         'ref_col_with_zero': [0.0, 2.0, 3.0, 4.0],
                         'ref_col_with_inf': [1.0, np.inf, 3.0, 4.0],
                         'ref_col_with_nan': [1.0, 2.0, np.nan, 4.0],
                         'target_col_with_inf': [np.inf, 4.0, 6.0, 8.0],
                         'target_col_with_nan': [1.0, 2.0, np.nan, 4.0]},
                         index=['r1', 'r2', 'r3', 'r4'])


@pytest.fixture(
        name="good_ref_col",
)
def good_ref_col():
    """Fixture for good reference column."""
    return 'reference_cond'


@pytest.fixture(
        name="good_target_col",
)
def good_target_col():
    """Fixture for good target column."""
    return 'target_cond'


@pytest.fixture(
        name="input_normalized_enzyme_activity",
)
def input_normalized_enzyme_activity():
    """Fixture normalized enzyme activity input."""
    return {'r1': 0.75, 'r2': 1.25, 'r3': 0.1, 'r4': 1.1}


@pytest.fixture(
        name="expected_dict_from_get_enzyme_bounds",
)
def expected_dict_from_get_enzyme_bounds():
    """Fixture for expected bounds from enzyme activity for output comparison."""
    return {'r1': 750.0, 'r2': 12.5, 'r3': 0.5, 'r4': 1100.0}
