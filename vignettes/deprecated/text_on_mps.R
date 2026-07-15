
# Designing custom management procedures

## Index ratio approaches

By default, the invertebrate MPs fish at current fishing effort levels without total allowable catch (TAC) control. You may however wish to invetigate options that explicity set catch limits. The simplest way to do this are empirical management procedures that use data directly and do not require an estimation / assessment step. These are increasingly applied in the management of a very wide range of fish species and their adoption is the principle objective of MSE in most fishery settings. 

The most commonly applied type of empirical  management procedure aims to fish at a stable exploitation rate by setting catch limits to a fixed multiple of a relative abundance index: 
  
  
  ## Index target approaches
  
  ## Index slope approaches
  
  In cases where you wish to maintain a trajectory in stock biomass at some level, an alternative to index target MPs are index slope MPs. Similarly to the index target approaches you must specify a target slope ('IS.targ') and two other arguments that control the calculation of the slope ('IS.yrs') and the responsiveness of the MP ('IS.fac'):
  
  ## Stock assessments 
  
  Management procedures can include estimation models that aim to mimick the real stock assessment processes. These are sometimes referred to as 'model-based MPs'. 

Even in situations where stock assessments are not possible for data, logistical or management reasons, it may be desirable to evaluate their potential performance (perhaps as a yardstick for current approaches). In MSE projections, all data are simulated including types that may not be currently available. That allows for the theoretical testing of alternative more data intensive options. 

The OpenMSE libraries come with a range of prespecified data-rich assessment MPs that include state-space surplus production (SP), state-space delay differential (DDSS) and statistical catch at age (SCA) models. A few derivatives are available that fish at constant FMSY ('_MSY'), constant 75% of FMSY ('_75MSY') and use a 40-10 HCR ('_4010', FMSY above 40% B0, zero fishing below 10% B0):
  
  ```{r args_MP_MC_datarich, eval=F}

```



## Short-cuts

Running model-based MPs can be computationally intensive. An alternative is to take the control points of an HCR not from an estimation model but from the operating model directly assuming some realistic level of observation (data going into the assessment) and estimation error (error in estimates of control points of an HCR in this context). Because there is no estimation, short-cuts can provide indicative results in a small fraction of the time as the full model-based MP. 

## Additional TAC control parameters

When using the specified MPs such as There are a number of other arguments that apply to all methods of catch advice (TAC.calc is set to 'Slope','Target' or 'Rate') that provide additional control over TAC changes. These include:
  
  * max TAC - that highest TAC level that is recommended by the MP
* minTAC - the lowest TAC level that is recommended by the MP
* TACdec - the highest fractional TAC decline that the MP can make (e.g. 0.2 means that TAC(y) must be higher than 0.8 x TAC(y-1))
* TACinc - the highest fraction TAC increase that the MP can make (e.g. 0.3 mean that TAC(y) must be lower than 1.3 * TAC(y-1))
* I.enp - Indices are processed via a polynomial smoother to try to remove noise from signal. The degree of smoothing is controlled by this parameter. ENP stands for effective number of parameters per datum. So higher ENPs are more closely fitting (less smooth) polynomials. The default of 0.25 means that for every 4 data points the polynomial adds another parameter for smoothing. 
* I_freq - the frequency that the indices will be observed in the forward projection. Let us assume there are three indices available a vector I_freq = c(1,3,0) would mean that there would be an observation every year for index 1, a new observation every three years for index 2, and no new data for index 3. 
* calib_yrs - the number of years used to define 'recent historical levels'. When arguments such as curI_2_target are used, they determine 'current' based on the calib_yrs number of years in the historical time period. So for example if the OMs were fitted to 2023 and calib_yrs was set to 3, the MP would use years 2021-2023 to determine 'recent' index and catch levels. 


## Other management levers

A default discarding rate and post-release mortality rate are also assumed. The shift to a new gear might be expected to change both of these. These are easily changed:
  
  
  ```{r Eff_DR_Fdisc, eval=F}

```

