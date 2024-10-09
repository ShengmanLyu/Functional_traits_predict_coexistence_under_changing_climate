# Functional traits predict coexistence under changing climate
Here are the data and R codes used in Lyu & Alexander (2024). Functional traits predict outcomes of current and novel competition under warmer climate. Global Change Biology. A full-text is availalbe on BioRxiv: https://doi.org/10.1101/2024.09.26.615168

Detailed information on population modelling and coexstence analyses is avaialbe:  Lyu, S. and J. M. Alexander (2023). "Compensatory responses of vital rates attenuate impacts of competition on population growth and promote coexistence." Ecology Letters 26(3): 437-447. https://doi.org/10.1111/ele.14167 R codes are availabel on GitHub at: 

### This repository includes data and R scripts used for population modelling and analyses in Lyu, S. & Alexander, J.M. (2023) Compensatory responses of vital rates attenuate impacts of competition on population growth and promote coexistence. Ecology Letters, 00, 1– 11. Available from: https://doi.org/10.1111/ele.14167

### There are four data files: 
#### Estimated_vital_rates.xlsx  
> This file includes the estimated vital rate parameters used for population modelling. Data and R codes used for estimating vital rate parameters can be found in Lyu, S. and Alexander, J. (2022). "Competition contributes to both warm and cool range edges." Nature Communications 13(1): 2502. https://doi.org/10.1038/s41467-022-30013-3

#### lambda_perturbed.RData  
> These are the projected population growth rates (lambda) of perturbed IPMs used for quantifying vital rate contributions. These can also be computed using R codes provided below. 

#### lambda_replaced.RData
> These are the projected population growth rates (lambda) of substituted IPMs when compensatory responses are present vs absent. These can also be computed using R codes  provided below. 

#### Coexistence_null.xlsx
> This is a template data sheet to store estimated niche overlap, relative fitness differences and the outcomes of competition.

### There are two R files:
#### Functions  
> R functions used for IPM, statistical analyses and data visualisation.

#### Modelling and analyses
> R scripts for statistical analyses and data visualisation (Figures 1-3 in the main text).
