"""Script to run the Eflux2 Algorithm."""

import cobra
import numpy as np
import pandas as pd
from optlang.symbolics import add
from eflux.utils import get_gpr_dict


def constrain_model_to_external_fluxes(
    model: cobra.Model, external_flux_scaling_factors: dict[str, float]
) -> cobra.Model:
    """Constrain model to external fluxes.

    inputs:
        model: cobra model
        external_flux_scaling_factors: dict (or dataframe column) of external fluxes for one strain/experimental condition scaled with respect to reference strain
    """
    # Get exchange reactions for external fluxes
    # FIXME: external_rxns being called here, but method is returning on model? if intentional, we probably don't need to be saving this to a variable.
    external_rxns = [
        rxn
        for rxn in model.reactions
        if rxn.id.startswith("EX_") and rxn.id[3:] in external_flux_scaling_factors
    ]

    # for rxn in model.reactions:
    #     if rxn.id in normalized_external_fluxes:
    #         rxn.lower_bound = normalized_external_fluxes[rxn.id]

    return model


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


def run_eflux(model: cobra.Model, tolerance: float = 1e-9) -> pd.DataFrame:
    """Run eflux: optimize the constrained cobra model by running FBA (perhaps with slack variables) to get a feasible solution.

    inputs:
        model: cobra model
        tolerance: tolerance for flux bounds
    outputs:
        eflux2_sol: dataframe of eflux2 solution
    """
    # Set model tolerance
    model.tolerance = tolerance

    # Constrain model to external fluxes
    constrained_model = model.copy()
    constrained_model = constrain_model_to_external_fluxes(
        constrained_model,
        {rxn.id: rxn.lower_bound for rxn in model.reactions if rxn.id.startswith("EX_")},
    )

    # Try running FBA to return a feasible solution
    try:
        # Constrain the model to enzyme activity

        # Succesful run with feasible solutions
        fba_sol = constrained_model.optimize(raise_error=True)
    # FIXME: Added the full path to OptimizationError, but seems like this still needs finished as e is not used.
    except cobra.exceptions.OptimizationError as e:
        # Use slack variables to overcome infeasibility
        slack_variables = {}

    # print("FBA status", fba_sol.status)
    # print("FBA solution", fba_sol.objective_value)

    # # Constrain the biomass to the optimal value
    # for r in constrained_model.reactions:
    #     if r.objective_coefficient:
    #         r.lower_bound = fba_sol.objective_value

    # # minimize the sum of squared flux values
    # constrained_model.objective = constrained_model.problem.Objective(
    #     add([r.flux_expression**2 for r in constrained_model.reactions]), direction="min"
    # )
    # eflux2_sol = constrained_model.optimize()
    # print("EFlux2 status", eflux2_sol.status)
    # print("EFlux2 solution", eflux2_sol.objective_value)
    # # return eflux2 solution
    # return eflux2_sol

    return eflux2_sol


def eflux2(
    constrained_model: cobra.Model,
    expression_scaling_factors: dict[cobra.Reaction, float],
    normalized_external_fluxes: dict[cobra.Reaction, float],
) -> pd.DataFrame:
    """Run eflux2.

    inputs:
        model: cobra model, already constrained to observed data for external fluxes & transcriptomics data

        normalized_enzyme_activity: dict (or dataframe column) of enzyme activity for one strain/experimental condition normalized to reference strain
        normalized_external_fluxes: dict (or dataframe column) of external fluxes for one strain/experimental condition normalized to reference strain
    outputs:
        eflux2_sol: dataframe of eflux2 solution
    """
    # Make a copy of the model
    # FIXME: Does this need to be constrained_model?
    eflux2_model = model.copy()

    # Get GPR dict containing isozymes (separated by 'or'),
    # where each isozyme has a set of subunits (separated by 'and')
    gpr_dict = get_gpr_dict(eflux2_model)

    # TODO: Move setting bounds to another function before running this one
    #
    # Set the bounds using the transcriptomics data
    for r in eflux2_model.reactions:
        if r.gene_reaction_rule:
            t = np.sum(
                [
                    np.min(
                        [
                            # FIXME: transcriptomics is not being passed into this method so this variable does not exist.
                            transcriptomics.loc[g]
                            if g in transcriptomics.index
                            else np.array([np.Inf])
                            for g in p
                        ]
                    )
                    for p in gpr_dict[r.id]
                ]
            )
            if r.lower_bound < 0.0:
                r.lower_bound = -t
            else:
                pass
            if r.upper_bound > 0.0:
                r.upper_bound = t
            else:
                pass
        else:
            if r.lower_bound <= -1000.0:
                r.lower_bound = -np.Inf
            if r.upper_bound >= 1000.0:
                r.upper_bound = np.Inf

    # solve FBA to calculate the maximum biomass
    eflux2_model.tolerance = 1e-9
    fba_sol = eflux2_model.optimize()
    print("FBA status", fba_sol.status)
    print("FBA solution", fba_sol.objective_value)
    # Constrain the biomass to the optimal value
    for r in eflux2_model.reactions:
        if r.objective_coefficient:
            r.lower_bound = fba_sol.objective_value
    # minimize the sum of squared flux values
    eflux2_model.objective = eflux2_model.problem.Objective(
        add([r.flux_expression**2 for r in eflux2_model.reactions]), direction="min"
    )
    eflux2_sol = eflux2_model.optimize()
    print("EFlux2 status", eflux2_sol.status)
    print("EFlux2 solution", eflux2_sol.objective_value)
    # return eflux2 solution
    return eflux2_sol
