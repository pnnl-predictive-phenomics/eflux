"""Utils module for eflux package."""

from typing import Tuple

import cobra
import numpy as np
import pandas as pd
from cobra import Gene, Reaction


def get_max_flux_bounds(
    model: cobra.Model, rxn_list: list[str], precision: int = 9
) -> Tuple[cobra.Model, pd.DataFrame]:
    """Get flux bounds from FVA to use in surrogate model of reference strain.

    Note: FVA = flux variability analysis
    inputs:
        model: cobra model
        rxn_list: list of reactions of interest, corresponding to reference strain selection criteria
        zero_threshold: magnitude threshold to identify and replace numerically zero flux values
    outputs:
        max_flux_bounds: max flux values to be used as a representative bounds of the reference strain.
    """
    # Run FVA to get (reasonably) tight bounds for all other reactions
    keep_rxn_list = [r.id for r in model.reactions if (r.id not in rxn_list)]
    flux_bounds = cobra.flux_analysis.flux_variability_analysis(
        model=model, reaction_list=keep_rxn_list, fraction_of_optimum=0.85, processes=8
    )
    max_flux_bounds = flux_bounds["maximum"].round(decimals=precision).to_dict()

    return max_flux_bounds


def get_gpr_dict(model: cobra.Model) -> dict[Reaction, set[frozenset[Gene]]]:
    """Gene reaction rule (GPR) for each reaction in the model.

    inputs:
        model: cobra model
    outputs:
        gpr_dict: dictionary of reactions to isozyme sets (corresponding genes from gene reaction rules)
    """
    # Parse GPR into a dict containing isozymes (separated by 'or')
    # Each isozyme has a set of subunits (separated by 'and')
    gpr_dict = {}
    for r in model.reactions:
        if r.gene_reaction_rule:
            isozymes = set()
            for isozyme in [isozyme.strip("() ") for isozyme in r.gene_reaction_rule.split(" or ")]:
                isozymes.add(frozenset(gene.strip("() ") for gene in isozyme.split(" and ")))
            gpr_dict[r] = isozymes

    return gpr_dict


def gene_expression_to_enzyme_activity(
    model: cobra.Model, gpr: dict[Reaction, set[frozenset[Gene]]], expression: dict[Gene, float]
) -> dict[Reaction, float]:
    """Map gene expression to enzyme activity inputs.

    inputs:
        model: cobra model
        gpr: dictionary of reactions (keys) to list of list of genes (values) for the correpsonding gene reaction rule.
        expression: dictionary of gene names (keys) to values from [likely] observed transcriptomics data.
    outputs:
        enzyme_activity: dictionary of reactions (keys) to corresponding isozyme activity from observed data (value).
    """
    enzyme_activity = {}
    for rxn in model.reactions:
        # Initialize enzyme_activity for this reaction to 0-value
        # Note: Converted to NaN-value IF this reaction doesn't have any genes in its gene reaction rule
        # Obvious example: Exchange/transport reactions don't have corresponding genes in their reaction rule, so that will take a NaN-value
        enzyme_activity[rxn] = 0.0

        if rxn in gpr:  # ensure rxn has a gene_reaction_rule defined
            for isozyme in gpr[rxn]:
                # Initialize isozyme_activity for this isozyme to infinity
                # Note: infinity-value is preserved IF this isozyme is not present in the observed transcriptomics data
                isozyme_activity = np.inf
                for gene in isozyme:
                    # ensure gene in the isozyme is included in observed data
                    if gene in expression:
                        isozyme_activity = np.min([isozyme_activity, expression[gene]])
                enzyme_activity[rxn] += isozyme_activity
        else:
            enzyme_activity[rxn] = np.nan

    return enzyme_activity


def convert_transcriptomics_to_enzyme_activity(
    transcriptomics_data: pd.DataFrame, model: cobra.Model
) -> pd.DataFrame:
    """Convert transcriptomics data to enzyme activity.

    inputs:
        transcriptomics_data: dataframe of transcriptomics data
        model: cobra model
        # gpr: dictionary of reaction ids (keys) to list of list of genes (values) for the correpsonding gene reaction rule
    outputs:
        enzyme_activity_df: dataframe of enzyme activity converted from transcriptomics data
    """
    # Initialize empty dataframe
    enzyme_activity_df = pd.DataFrame()

    # Get gene production rules
    gpr = get_gpr_dict(model)

    # Loop through each strain to convert each column of transcriptomics data
    for this_strain in transcriptomics_data.columns:
        # Create dict of genes and corresponding float values using trancsciptomics data
        expression_dict = {
            g: transcriptomics_data.loc[g][this_strain] for g in transcriptomics_data.index
        }

        # Run the gene expression to enzyme activity converter for this_strain
        enzyme_activity_dict = gene_expression_to_enzyme_activity(model, gpr, expression_dict)

        # Initialize empty dataframe
        if this_strain == transcriptomics_data.columns[0]:
            # Use enzyme_activity_dict keys as the index
            enzyme_activity_df = enzyme_activity_df.reindex(enzyme_activity_dict.keys())
            # Add reaction ID column
            enzyme_activity_df["Reaction_ID"] = [k.id for k in enzyme_activity_dict]

        # Add enzymze_activity to dataframe
        enzyme_activity_df[this_strain] = enzyme_activity_dict

    if enzyme_activity_df.empty:
        return pd.DataFrame()

    return enzyme_activity_df.set_index("Reaction_ID")
