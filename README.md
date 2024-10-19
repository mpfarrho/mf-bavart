## Description
"Nowcasting in a Pandemic using Non-Parametric Mixed Frequency VARs" by Huber, F., Koop, G., Onorante, L., Pfarrhofer, M., and J. Schreiner, _Journal of Econometrics_, **232**(1), 2023, 52-69.
- [Publication](https://doi.org/10.1016/j.jeconom.2020.11.006)
- [Working paper](https://arxiv.org/abs/2008.12706)

### Code files
These files create a function mfbavart(...) to estimate the MF-BAVART model.
- mfbavart_func.R contains the main function
- aux_func.R collects several auxiliary functions
- example.R contains an example code for using the function.

In addition to the baseline model in the paper, the code also includes an option to introduce stochastic volatility (SV) in the error terms. Several parts of the original code used for the paper in directory "replication" have been replaced to improve computational efficiency.

This code comes without technical support of any kind. The code is free to use, provided that the paper is cited properly.

Some codes and helper functions are taken or adapted from the "mfbvar" package. Thanks to Vincent Dorie (mtn. of "dbarts") and Sebastian Ankargren (mtn. of "mfbvar") for technical support regarding their excellent packages.

### Inputs for mfbavart(...):
- data            a list that contains ts-objects of different frequencies in its 
                  M (number of endogenous variables) slots, such that high-frequency (monthly)
                  series are ordered first and followed by low-frequency series (quarterly)
- itr             intertermporal restriction ("lvl" or "grw") of length 
                  corresponding to number of low frequency series
- p               numeric, lag-length of the VAR (minimum of 5 if itr=="grw", and 3 if itr=="lvl")
- fhorz           numeric, forecast horizon in months (3 per quarter)
- cons            TRUE/FALSE, whether a constant should be included
- exact           TRUE/FALSE, whether BART fit is stored in output values or 
                  filtered data based on approximation
- sv              TRUE/FALSE, whether structural errors feature stochastic volatility
- var.thrsh       numeric, threshold for resampling coefficients by draw (for sampler stability)
- max.count.var   numeric, maximum number of resampling steps
- cgm.level       numeric \in (0,1), \alpha in the paper (probability of terminal node)
- cgm.exp         numeric > 0, \beta in the paper (probability of terminal node)
- sd.mu           numeric, \gamma in the paper
- num.trees       numeric, number of trees for BART, S in the paper
- prior.sig       numeric of length 2, [1] nu_j, [2] v in the paper, 
- nburn           numeric, number of burnins
- nsave           numeric, number of draws for posterior/predictive inference
- thinfac         numeric, thinning factor
- quiet           TRUE/FALSE, whether progress bar is indicated during sampling

Function returns:
- Y              [nsave,T,M] array that contains the latent states (see also option "exact")
- fcst           [nsave,fhorz,M] array that contains forecasts
- Yq             [nsave,T+fhorz,M] array that contains aggregated filtered and forecasted series
