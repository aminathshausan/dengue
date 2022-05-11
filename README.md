# Minimizing severity of dengue serotype 1 infection by transmissible interfering particles

This repository makes available the source code for the work: "Minimizing severity of dengue serotype 1 infection by transmissible interfering particles", available online as an open access preprint on bioRxive: doi: https://doi.org/10.1101/2020.04.21.052936.

# Prerequisites

To run the codes you will require R(>=3.30) and latest version of RStudio.

# Directory structure
```
├── project
│   ├── data
│   ├── images
│   ├── results
│   ├── src
│   │   ├── plot_raw_data_and_sensivity.R
│   │   ├── run_models_sim.R
│   │   ├── analyze_fits.R
│   │   ├── TIPs_therapy.R
|   |   ├── model_a_sim.stan
|   |   ├── model_b_sim.stan
```
- The `data` folder contains all preprocessed datasets reqruired to run the codes in the `src` folder. 
- The `images` folder contains plots generated from the codes in `src` folder.
- The `results` folder contains all results saved (as .RData format) after runing the codes in the `src` folder.
- The `src` folder contains all source codes required to produce the results and images.

# How to run scripts

1. Run the script `plot_raw_data_and_sensivity.R` to generate plots for `Fig 1` and `Fig 2` shown in the `images` folder
2. Run the script `run_models_sim.R` to fit the Bayesian hierarchical model `model_a_sim.stan` or `model_b_sim.stan`. These models are written using the Stan programing language.
3. Run the scrip `analyze_fits.R` to analyze results of the fits and to generate plots for Figures 3-4 S1, S2 and S5 shown in the `images` folder.
4. Run the scrip `TIPs_therapy.R` to simulate results for experiments 1, 2 and the individual effects.  This script is also used to generate plots for Figures 5-7, S3, S4 and S6-S16 shown in the `images` folder.
