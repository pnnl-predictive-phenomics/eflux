# E-Flux (combining Expression Data with Fluxes)

[![Actions status](https://github.com/pnnl-predictive-phenomics/eflux/workflows/CI/badge.svg)](https://github.com/pnnl-predictive-phenomics/eflux/actions)
[![Rye](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/rye/main/artwork/badge.json)](https://rye-up.com)
[![License](https://img.shields.io/badge/License-BSD_2--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

Package to run E-Flux on `cobra` models with provided transcriptomics expression data.
# eflux

A repository for hosting the code for E-flux.

```mermaid
flowchart TD
 subgraph Legend[<h2> Legend </h2>]
 direction LR
        UserData["User Provided Data"]
        CompData["Computed Data"]
        Action(("Action/Process"))
  end

 subgraph ObservedData[<h2> Observed Data and Pre-processing </h2>]
        Data["Observed Data"]
        TranscData["Transcriptomics"]
        ExMetab["External Metabolites"]
        RxnConstr["Reactions of Interest & Objectives"]
        CalcExFlux["Calculated External Fluxes"]
        EnzAct["Enzyme Activity"]
  end

 subgraph CobraModel[<h2> Cobra Model and Pre-processing </h2>]
        FVA(("FVA"))
        OrigCobra["Cobra Model"]
        FluxBounds["Flux Bounds"]
        RefCobra["Data-Conditional Cobra Model"]
        DelZeroFlux(("Delete 0-Flux Reactions"))
  end

  subgraph E-flux[<h2> E-flux Algorithm </h2>]
        Eflux(("Eflux3"))
        CondCompFluxes["Conditional Computed Fluxes"]
  end

    %% Data Links
    Data --> TranscData & ExMetab & RxnConstr 
    ExMetab --> CalcExFlux
    TranscData -- Convert --> EnzAct

    %% Cobra Model Links
    OrigCobra --> FVA 
    FVA --> FluxBounds
    FluxBounds --> RefCobra
    RefCobra --> DelZeroFlux & Eflux 
    DelZeroFlux --> RefCobra
    Eflux --> CondCompFluxes

    %% Inter-subgraph Links
    RxnConstr --> FVA
    CalcExFlux & EnzAct --> Eflux

    %% Color styles for each node
    %% Raw Data - Blue
    style UserData stroke:#2962FF,fill:#2962FF,color:#FFFFFF
    style TranscData stroke:#2962FF,fill:#2962FF,color:#FFFFFF
    style ExMetab stroke:#2962FF,fill:#2962FF,color:#FFFFFF
    style RxnConstr stroke:#2962FF,fill:#2962FF,color:#FFFFFF
    style Data stroke:#2962FF,fill:#2962FF,color:#FFFFFF
    style OrigCobra stroke:#2962FF,fill:#2962FF,color:#FFFFFF

    %% Computed Data- Green
    style CompData fill:#00C853,color:#000000
    style CalcExFlux fill:#00C853,color:#000000
    style EnzAct fill:#00C853,color:#000000
    style FVA stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,fill:#FFD600,color:#000000
    style FluxBounds fill:#00C853,color:#000000
    style RefCobra fill:#00C853,color:#000000
    style CondCompFluxes fill:#00C853,color:#000000

    %% Action/Process - Yellow
    style Action stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,fill:#FFD600,color:#000000
    style DelZeroFlux stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,fill:#FFD600,color:#000000
    style Eflux fill:#FFD600,stroke:#FFD600,stroke-width:4px,stroke-dasharray: 0,color:#000000
```

## Getting Started 🏃

To get started run the following:

```python
from eflux import eflux2

eflux2(cobra_model, transcriptomics)
```

## Installation 🪛

The most recent code and data can be installed directly from GitHub with:

```shell
pip install git+https://github.com/pnnl-predictive-phenomics/eflux.git
```

## License 📄

The code in this package is licensed under the BSD-2 License.


## Contributing 👋
To contribute to this package, please reference [CONTRIBUTING.md](CONTRIBUTING.md)
