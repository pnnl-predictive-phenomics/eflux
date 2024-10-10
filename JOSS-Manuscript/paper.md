---
title: 'E-flux: A Python package for calculating metabolic fluxes conditioned on gene expression'
tags:
  - Python
  - 
  - 
  - 
  - 
authors:
  - name: Shant M. Mahserejian
    orcid: 0000-0001-5360-6039
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Andrew McNaughton
    orcid: 0000-0002-4146-7921
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: Jeremy D. Zucker
    orcid: 0000-0002-7276-9009
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
affiliations:
 - name: Pacific Northwest National Laboratory
   index: 1
   ror: 
#  - name: Institution Name, Country
#    index: 2
#  - name: Independent Researcher, Country
#    index: 3
date: 30 September 2024
bibliography: paper.bib

# # Optional fields if submitting to a AAS journal too, see this blog post:
# # https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Here's a citation to the first publication for [E-flux](https://doi.org/10.1371/journal.pcbi.1000489) [@eflux1].

Eflux1 quote: "Our method, which we call E-Flux (as a combination of flux and expression), extends the technique of Flux Balance Analysis by modeling maximum flux constraints as a function of measured gene expression. In contrast to previous methods for metabolically interpreting gene expression data, E-Flux utilizes a model of the underlying metabolic network to directly predict changes in metabolic flux capacity.... E-Flux thus provides a promising new approach for algorithmically predicting metabolic state from gene expression data. "



# Introduction 

Flux Balance Analysis (FBA) is a widely used computational method in systems biology for predicting the flow of metabolites through metabolic networks based solely on stoichiometric constraints. Despite its utility, traditional FBA does not account for dynamic regulatory factors such as enzyme expression levels, which play a crucial role in modulating metabolic fluxes. This limitation can lead to less accurate predictions of metabolic behaviors under varying genetic and environmental conditions.

The **E-flux** software addresses this limitation by integrating enzyme expression data directly into FBA to enhance the accuracy of predicted metabolic fluxes. The core idea of E-flux is simple: enzyme concentrations constrain the maximum allowable fluxes through metabolic pathways, and changes in enzyme levels are linearly correlated with changes in these maximum fluxes. However, the implementation of this concept requires sophisticated techniques to translate relative enzyme expression measurements into absolute flux constraints.

To address this, we have developed three variants of the E-flux software, each designed to handle different computational and biological challenges:

**E-flux 1:** Utilizes an L1 norm to minimize the sum of the fluxes, promoting sparsity in the flux distribution.
**E-flux 2:** Employs an L2 norm to minimize the sum of the fluxes, promoting a more evenly distributed flux across pathways.
**E-flux 3:** Introduces slack variables to ensure flux constraints do not cause infeasibilities, thus providing a more robust solution when dealing with inconsistent or noisy enzyme data.



# Statement of need

Previous implementations of the Eflux algorithm aimed to extending flux balance analysis (FBA) to model metabolic fluxes conditioned on gene expression [@eflux1], [@eflux2]. This allowed for predicting fluxes based on data for different experimental conditions and gene expression profiles [@eflux1_usage]. However, some conditions can strain the model optimization process used to calculate fluxes, leading to infeasible solutions. To this end, we developed Eflux3, a new version of the algorithm that incorporates slack variables to ease the constraints in the optimization process and thus improve outcomes by allowing for a set of agreeable flux values.

Here we introduce Eflux as a Python package that allows for the calculation of fluxes conditioned on gene expression data. The package is meant to house the latest implementation of Eflux (Eflux3), as well as legacy implementations for completness and to allow for reproducing results.


# Software Description

## Overview
- Name: Introduce Eflux and its main purpose.
- Legacy Components: Describe the inclusion of Eflux1 and Eflux2 scripts.
- Novel Contribution: Introduce Eflux3 and its unique feature—incorporation of slack variables to improve FBA outcomes.

## Key Features
- User-Friendly Design: Highlight documentation, unit tests, and package structure.
- Solution Robustness: Explain how slack variables enhance the solver by addressing infeasible solutions.

## Implementation
- Technologies and Frameworks: Mention the use of Python and any other relevant technologies.
- Code Structure: Provide an overview of the package structure and key modules.

## Mathematics
blah blah
<!-- 
Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations: 

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text. -->


# Installation and Usage
- Installation Instructions: Guide on how to install the package.
- Basic Usage: Provide basic examples or command-line instructions to demonstrate the typical usage of Eflux.
- Advanced Features: Briefly touch on more advanced functionalities and configurations.


# Software Validation
- Unit Tests: Mention the incorporation of unit tests to ensure code quality and reliability.
- Documentation: Describe the comprehensive documentation available for users.


# Impact
- (Note: Keep focus on the software)
- Expected Use Cases: Briefly mention scenarios where Eflux can be particularly valuable to the research community.


# Acknowledgments
- Contributors: Thank those who contributed to the development of the package.
- Funding: This research was supported by Predictive Phenomics Initiative from October 1, 2022-September 30, 2024, under the Laboratory Directed Research and Development (LDRD) Program at Pacific Northwest National Laboratory (PNNL).  PNNL is a multi-program national laboratory operated for the U.S. Department of Energy (DOE) by Battelle Memorial Institute under Contract No. DE-AC05-76RL01830.  Funding after October 1, 2024 was supported by  the DOE Biosystems Design Program under award number DE-SC0023091.



# References

 

<!-- Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)" -->

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }


