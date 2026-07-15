# mahi MP
# Feb 22nd 2026

# Cassidy Notes:
# Here is a link to the current regulations: https://safmc.net/species/dolphin/, which includes tabs for rec vs. commercial regulations by state.
# Rec regulations: bag limit = 10 fish; vessel limit = 54 fish (effective May 2022); 20"FL min size in SC, GA, and FL
# Com regulations: Commercial ACL = 1719953 lbs Whole Weight
# John H would be a great person to reach out to with management questions.

#' The Mahi Management Procedure (OpenMSE v2.0)
#'
#' A function for emulating status quo mahi management and allowing for a wide range of combined management controls. Fleets (n=8) are: (1) USCom, (2) RecN, (3) RecS, (4) HireN, (5) HireS, (6) Intl, (7) Disc, (8) UnRep. Areas (n = 5) are:  (1) "CAR+FLK", (2) "NCA", (3) "NCFL", (4) "NED", (5) "NNC+VBM"
#'
#' @param Data An object of openMSE class 'Data'
#' @param TAC Annual total allowable catch. A single value (all fleets, all areas), a vector length nFleet or a matrix nFleet x nArea. Positive real number.
#' @param rel_TAC A single value that is current TAC relative to the most recent  historical year. Positive real number. Can be a single number, a vector nfleet long or a matrix nfleet by narea long.
#' @param TAE Annual total allowable effort. Single value (all fleets, all areas) or a vector length nFleet. Positive real number.
#' @param rel_TAE Relative effort (an improper fraction of effort in the last historical year). Positive real number. Single number (all fleets/areas), vector of length nfleets (8), vector of length nareas, or matrix with nfleets rows and nareas columns (8 x 7).
#' @param SLmin Minimum size limit (mm). Positive real number. Can be a single number (all fleets) or a vector nfleets long (8).
#' @param SLmax Maximum size limit (mm). Positive real number larger than Smin (if specified). Can be a single number (all fleets) or a vector nfleets long (8).
#' @param TL Trip limit (fish retained per trip). Positive integer. Can be a single number (all fleets) or a vector nfleets long (8). Note that predictive models are only available for RecS and HireS in regions CAR+FLK and NCFL.
#' @param PRM Post Release Mortality (discard mortality). List nfleet long of either a Fraction (all length classes) or a vector nlengths long. Default from US hook and line estimates of Rudershausen et al 2019 (DOI: 10.1002/nafm.10348). Can be a single number or a vector of length nareas (7).
#' @param HCR_use Boolean vector nHCR long. Should HCR(s) be used?  If False, all other HCR arguments will be ignored.
#' @param HCR_ICP List (nHCR long) of 2 - position vectors. The Index control points. These are a pair of x-axis control points for a hockey-stick control rule. Phrased in terms of the levels of the index (HCR_input) over the HCR_calib years.
#' @param HCR_LCP  List (nHCR long) of 2 - position vectors. Management Lever control points. These are a pair of y-axis control points for a hockey-stick control rule. Phrased in terms of the levels of the Lever (HCR_lever) over the last historical year (see objects rTAC & rEffort).
#' @param HCR_fleet The fleet for which the HCR is used. Vector of positive integers nHCR long. If NA this sets advice for all fleets combined based on the HCR_input specified.
#' @param HCR_sens Sensitivity of the response. A vector of positive real numbers nHCR long. A value of 1 is proportional.
#' @param HCR_lever The output variable. A vector of character strings nHCR long. E.g., "TAC", "TAE".
#' @param HCR_up_max The maximum rate of advice increase. Positive imperfect fraction. E.g., 0.2 is a maximum management increase of +20 per cent.
#' @param HCR_down_max The maximum rate of advice decrease. Positive fraction. E.g. 0.8 is a maximum management decreate of 80 per cent.
#' @param HCR_input The indices used in the harvest control rule. Vector of integers, nHCR long. E.g., HCR_input = c(2, 1) uses index 2 for the first HCR, index 1 for the second HCR
#' @param HCR_smooth The effective number of parameters of a polynomial smoother applied to the input indices. Defaults to NA - raw data. A value of 0.1 means that there will be length(Index)*0.1 number of smoothing parameters. Hence, higher values are less smoothing. This parameterization keeps smoothing consistent as the length of the time series increases in projections.
#' @param Index_names Vector of character vectors nIndex long. The indices to be generated. See object 'Indices'. The indices simulated. e.g. c("PLL", "PR_Tourn")
#' @param Index_cv  Vector of positive real numbers nIndex long. User override for the coefficient of variation of the index (CV). NA means stat properties are estimated from model fit
#' @param Index_ac Vector of fractions nIndex long. User override for the The lag-1 autocorrelation of the index. (Generally positive). NA means stat properties are estimated from model fit
#' @param Index_lag Vector of integers, nIndex long. The number of years of data lag. 1 is a single year.
#' @param Index_seed Real number vector nI long. A random number seed for generating index observation errors nIndex long
#' @param Index_beta Boolean vector nI long. Should the beta parameter be estimated in the index - biomass relationship I = qB^beta ? Strongly recommend NO!
#' @param TL_empirical Should the empirical distribution of catch rates be used to calculate release rate according to trip limits?
#' @param plot Should plots be produced that explain MP calculations (e.g. HCRs)? Boolean.
#' @param debugfile Character string, a file and directory to write an rds list of internal data e.g. "C:/temp/debug.rds". If 'NA', nothing is written
#' @param debugts Positive integer greater than number of historical years, the year that debug data are written.
#' @param debugDatafile Character string, a file and directory to write the Data object to e.g. "C:/temp/debugData.rds". If 'NA', nothing is written
#' @param docheck Conduct an internal test of consistency among MP arguments
#' @examples
#' mahiMP(Example_Data)
#' @author T. Carruthers
#' @export
mahiMP = function(Data,
                  TAC = NA,
                  rel_TAC = 1,
                  TAE = NA,
                  rel_TAE = 1,
                  SLmin = c(0,    0,    500,  0,    500,  0,    0,    0),
                  SLmax = c(1E4,  1E4,  1E4,  1E4,  1E4,  1E4,  1E4,  1E4),
                  TL =    c(NA,   54,   54,   54,   54,   NA,   NA,   NA),
                  PRM =  list(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25),
                  HCR_use = F,
                  HCR_ICP = NA,
                  HCR_LCP = NA,
                  HCR_fleet = list(NA),
                  HCR_sens = 1,
                  HCR_lever = "TAC",
                  HCR_up_max = 0.5,
                  HCR_down_max = 0.5,
                  HCR_input = 1,
                  HCR_smooth = 0.3,
                  Index_names = NA,
                  Index_cv = NA,
                  Index_ac = NA,
                  Index_lag = 1,
                  Index_seed = 1,
                  Index_beta = F,
                  TL_empirical = T,
                  plot = F,
                  onlyHCR = F,
                  debugfile = NA,
                  debugts = 152,
                  debugDatafile = NA,
                  docheck = F,
                  ...){

  # ARGS:
  #   Data = Example_Data_2; TAC = NA; rel_TAC = NA; TAE = NA; rel_TAE = 1;  SLmin = c(0,    0,    500,  0,    500,  0,    0,    0);  SLmax = c(1E4,  1E4,  1E4,  1E4,  1E4,  1E4,  1E4,  1E4)
  #   TL =    c(NA,   54,   54,   54,   54,   NA,   NA,   NA);  PRM =  list(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25)
  #   HCR_use = F; Index_names = "PLL";  Index_cv = NA;  Index_ac = NA;  Index_lag = 1; Index_seed = 1; Index_beta = F; TL_empirical = T; plot =T
  #   debugfile = NA; debugts = 158; debugDatafile = NA; docheck = F


  mahiMPcheck(docheck, Data, TAC, TAE, SLmin, SLmax, TL)             # Check dimensions and specs [0]

  Sim_Ind = gen_ind_2(Data, Index_names, Index_cv, Index_ac,         # Generate abundance indices [2]
                      Index_lag, Index_seed, Index_beta, plot=plot)

  season = get_season(Data)                                          # Get season (this time step, not last with data) [14]

  ad = new('advice')                                                 # Make new Advice object
  ad = do_TAC(Data, TAC, rel_TAC, season, ad)                        # Add TAC info if applicable [55]
  ad = do_TAE(Data, TAE, rel_TAE, season, ad)                        # Add TAE info if applicable [29]
  ad = do_PRM(Data, PRM, ad)                                         # Adds the post release mortality [1]
  ad = do_size(Data, SLmin, SLmax, ad)                               # Adds any specified size limits [2]
  ad = do_TL(Data, TL,  plot, TL_empirical, ad)                      # Apply any bag limits where there are data BL_info [13]
  ad = do_HCR(Data, HCR_use, HCR_ICP, HCR_LCP, HCR_fleet, HCR_sens,  # Harvest control rule comes last as it adds dynamic rules [762]
              HCR_lever, HCR_up_max, HCR_down_max, HCR_input,
              HCR_smooth, Sim_Ind, plot = plot, onlyHCR = onlyHCR, ad)
  # ad = do_interim()                                                # Interim (post conditioning, pre MP) data

  if(plot) plot_mahi_ad(ad)                                          # Plots internal calculations (smoothing, TL, HCR)

  writedebug(debugfile, debugts, Data, Sim_Ind, debugDatafile)       # Internal debugging data if debugfile !=NA [0]

  ad
}

class(mahiMP) = "mp"



