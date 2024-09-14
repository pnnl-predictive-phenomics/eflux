"""Script to run the Eflux2 Algorithm."""

from typing import Optional

import cobra
import numpy as np
import pandas as pd

from eflux.utils import get_max_flux_bounds, load_model_from_path


def add_slack_variables_to_model(
    model: cobra.Model, upper_bounds: dict[str, float], slack_weight: float = 1000
) -> cobra.Model:
    """Add slack variables to model.

    inputs:
        model: cobra model with objective already defined (TBD, maybe as optional.....and upper bounds defined by FVA)
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
        this_slack_var = relaxed_model.problem.Variable("SLACK_" + r_id, lb=0)
        slack_vars.append(this_slack_var)
        # Add a constraint between the reaction flux and the slack variable using the upper bound
        constraint = relaxed_model.problem.Constraint(
            this_rxn.flux_expression - this_slack_var, lb=0, ub=bound
        )
        relaxed_model.add_cons_vars(constraint)

    # Define a new combined objective
    combined_objective = relaxed_model.problem.Objective(
        relaxed_model.objective.expression - slack_weight * sum(slack_vars), direction="max"
    )

    # Set the combined objective as the objective of the model
    relaxed_model.objective = combined_objective

    return relaxed_model


def get_normalized_condition(
    df: pd.DataFrame, *, ref_col: str, target_col: str
) -> dict[str, float]:
    """Create dictionary of normalized/relative scale factors for a target condition w.r.t. a reference condition.

    inputs:
        df: observed data with reaction ids as rownames, and column names of experimental conditions that include the reference and target
        ref_col: column name for the reference condition
        target_col: column name for the target condition
    ouptut:
        norm_cond_dict: a dictionary of reaction ids (keys) mapped to scaling factors (values)
    """
    # Check zero or inf or nan entries in ref_col
    if np.any(df[ref_col] == 0) or np.any(np.isinf(df[ref_col])) or np.any(np.isnan(df[ref_col])):
        raise ValueError("Reference condition contains zero or inf or nan entries")

    # Check for inf or nan entries in target_col
    if np.any(np.isinf(df[target_col])) or np.any(np.isnan(df[target_col])):
        raise ValueError("Target condition contains inf or nan entries")

    # Calculate relative scale factors
    norm_cond_dict = df[target_col].divide(df[ref_col]).to_dict()

    return norm_cond_dict


def get_condition_specific_upper_bounds(
    fva_upper_bounds: dict[str, float], scaling_factors: dict
) -> dict[str, float]:
    """Get upper bounds for one experimental condition/strain only.

    inputs:
        fva_upper_bounds: dictionary of reaction ids (keys) and upper bounds from FVA (values)
        scaling_factors: dict of of reaction id (keys) and scaling_factors (values) obtained by
                        normalizing observed data for one strain/experimental condition with
                        respect to a reference condition
    outputs:
        dict of reaction id (keys) and upper bound on model reaction fluxes (values) for corresponding to one strain/experimental condition
    """
    return {r: b * scaling_factors[r] for r, b in fva_upper_bounds.items() if r in scaling_factors}


def run_condition_specific_eflux(
    model: cobra.Model,
    growth_rxn_id: str,
    product_rxn_id: str,
    external_fluxes: pd.DataFrame,
    enzyme_activity: pd.DataFrame,
    ref_cond: str,
    target_cond: str,
) -> dict[str, float]:  # or -> pd.DataFrame
    """Run eflux for one strain/experimental condition.

    inputs:
        model: cobra model with objective already defined
        growth_rxn_id: reaction id for growth
        product_rxn_id: reaction id for product
        external_fluxes: dataframe of external fluxes
        enzyme_activity: dataframe of enzyme activity
        ref_cond: reference condition (column of both external_fluxes and enzyme_activity)
        target_cond: target condition (column of both external_fluxes and enzyme_activity)
    outputs:
        fluxes: dictionary of reaction ids (keys) and flux values (values)
    """
    return {}

    # def adjust_reaction_bounds(model: cobra.Model, flux_bounds: pd.DataFrame) -> cobra.Model:
    #     """Adjust reaction bounds by FVA."""
    #     new_model = model.copy()

    #     for r in flux_bounds.index:
    #         new_model.reactions.get_by_id(r).lower_bound = 0
    #         new_model.reactions.get_by_id(r).upper_bound = flux_bounds.loc[r, "maximum"]

    #     return new_model

    # def prep_model_for_eflux(model: cobra.Model) -> cobra.Model:
    """Check model feasiblity and adjust bounds using FVA."""
    # Check if model is feasible

    # Check model growth rate

    # Check model production rate

    # Set model tolerance
    # model.tolerance = tolerance

    # Run FVA to get flux bounds for all reactions,  external reactions


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


def eflux3(
    model_path: str, data_path: str, *, reference_col: str, objective: Optional[str] = None
) -> pd.DataFrame:
    """Run the Eflux3 workflow.

    Parameters
    ----------
        model_path : str
            The path to the CobraPy compliant model for running eflux on.
        data_path : str
            The path to the transcriptomics data to condition the fluxes on.
        reference_col: str
            The column name of the reference column.

    Returns
    -------
        df : pd.DataFrame

    """
    enzyme_activity = pd.read_csv(data_path, index_col="Reaction_ID")

    model = load_model_from_path(model_path)
    # model.reactions.BIOMASS__1.lower_bound = 0.01
    # model.reactions.EX_sucr_e.lower_bound = 0.0

    if objective is not None:
        model.objective = objective

    fva_upper_bounds: dict[str, float] = get_max_flux_bounds(
        model,
        rxn_list=[],
    )  # ['EX_sucr_e', 'EX_photon650_e', 'EX_photon690_e', 'EX_co2_e', 'BIOMASS__1'])

    condition_specific_upper_bounds = {
        condition: get_condition_specific_upper_bounds(
            fva_upper_bounds, enzyme_activity[condition].dropna().to_dict()
        )
        for condition in enzyme_activity
    }
    condition_specific_models = {
        condition: add_slack_variables_to_model(model, upper_bounds)
        for condition, upper_bounds in condition_specific_upper_bounds.items()
    }

    condition_specific_fluxes = pd.DataFrame(
        {
            condition: condition_specific_model.optimize().fluxes
            for condition, condition_specific_model in condition_specific_models.items()
        }
    )

    fluxes = condition_specific_fluxes
    # fluxes.to_csv('../data/circadian_experiments/processed_data/enzyme_constrained_fluxes.csv')
    fluxes.to_csv(
        "../data/circadian_experiments/processed_data/enzyme_constrained_fluxes_nocofactorsinmodel.csv"
    )
