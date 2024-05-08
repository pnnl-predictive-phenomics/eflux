"""Script to run the Eflux2 Algorithm."""

import cobra
import numpy as np
import pandas as pd
from optlang.symbolics import add

from eflux.utils import get_gpr_dict


def add_slack_variables_to_model(model: cobra.Model, upper_bounds: dict[str, float], slack_weight: float = 1000) -> cobra.Model:
    """Add slack variables to model.

    inputs:
        model: cobra model with objective already defined
        upper_bounds: dict (or dataframe column) of reaction id keys and upper bound values for fluxes corresponding to one 
                      strain/experimental condition (e.g. scaled/normalized enzyme activity or external fluxes)
        slack_weight: weight of slack variables relative to model.objective
    outputs:
        model: cobra model constrained using upper bounds, but relaxed using slack variables
    """
    # Check for correct input types
    if model is None:
        raise TypeError("model cannot be None")
    if upper_bounds is None:
        raise TypeError("upper_bounds cannot be None")

    # Copy model to prevent overwriting
    relaxed_model = model.copy()

    # Initialize list of slack variables
    slack_vars = []

    # Add slack variables
    for r_id, bound in upper_bounds.items():
        # Create slack variable for each reaction
        this_rxn = relaxed_model.reactions.get_by_id(r_id)
        this_slack_var = relaxed_model.problem.Variable('SLACK_' + r_id, lb=0)
        slack_vars.append(this_slack_var)
        # Add a constraint between the reaction flux and the slack variable using the upper bound
        constraint = relaxed_model.problem.Constraint(this_rxn.flux_expression - this_slack_var, lb=0, ub=bound)
        relaxed_model.add_cons_vars(constraint)

    # Define a new combined objective
    combined_objective = relaxed_model.problem.Objective(
        relaxed_model.objective.expression - slack_weight * sum(slack_vars),
        direction='max'
    )

    # Set the combined objective as the objective of the model
    relaxed_model.objective = combined_objective

    return relaxed_model


def get_enzyme_bounds(model: cobra.Model, norm_enzyme_activity: dict) -> dict[str, float]:
    """Get enzyme bounds.

    inputs:
        model: cobra model with reaction bounds adjusted by FVA
        norm_enzyme_activity: dict (or dataframe column) of enzyme activity for one strain/experimental condition normalized with respect to reference strain
    outputs:
        enzyme_bounds: dict of reaction id keys and upper bound values for fluxes corresponding to one strain/experimental condition
    """
    bounds_dict = {}



    return {}


def get_normalized_condition(df:pd.DataFrame, ref_col: str, taget_col: str) -> dict[str, float]:
    """"""
    # Check zero or inf entries in ref_col

    # Check for inf entries in target_col

    return df[taget_col].divide(df[ref_col]).to_dict()


# Main function expected flow:
# 1) inputs: cobra model, reference condition, observed external fluxes (can be metabolomics), observed enzyme activity (can be transcriptomics)
# 2) run FVA
# 3) Adjust cobra model bounds using FVA bounds
# 4) Optional: convert metabolomics to external fluxes
# 5) Optional: convert transcriptomics to enzyme activity
# 6) Normalize external fluxes wrt reference condition
# 7) Normalize enzyme activity fluxes wrt reference condition
# 8) Get bounds for all reactions with corresponding observed data (from both enzyme activity and external fluxes)
# 9) Run model constraints while using slack variables (add_slack_variables_to_model)
# 10) Run FBA and output fluxes
