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
