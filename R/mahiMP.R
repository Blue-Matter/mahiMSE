# mahi MP
# July 25th 2025
# Currently in development. Delivery date August 2025.
# Notes: needs to be updated to FSRA (seasonal F matrices)

# Hi Tom --

# Cassidy Notes:
# Here is a link to the current regulations: https://safmc.net/species/dolphin/, which includes tabs for rec vs. commercial regulations by state.
# Rec regulations: bag limit = 10 fish; vessel limit = 54 fish (effective May 2022); 20"FL min size in SC, GA, and FL
# Com regulations: Commercial ACL = 1719953 lbs Whole Weight
# John H would be a great person to reach out to with management questions.

#' The Mahi Management Procedure
#'
#' A function for emulating status quo mahi management and allowing for a wide range of combined management controls. Fleets (n=8) are: (1) USCom, (2) RecN, (3) RecS, (4) HireN, (5) HireS, (6) Intl, (7) Disc, (8) UnRep. Areas (n = 7) are:  (1) NED, (2) NE, (3) NC, (4) SE, (5) SFL, (6) SAR, (7) CAR.
#'
#' @param x Simulation number
#' @param DataList A list object [[nstocks]][[nareas/fleets]] of openMSE objects of class 'Data'
#' @param Effort Relative effort (an improper fraction of effort in the last historical year). Positive real number. Single number (all fleets/areas), vector of length nfleets (8), vector of length nareas, or matrix with nfleets rows and nareas columns (8 x 7).
#' @param TAC Annual TAC (kg). Positive real number. Single number (all fleets/areas), vector of length nfleets (8), vector of length nareas (7), matrix with nfleets rows and nareas columns (8 x 7).
#' @param Smin Minimum size limit (mm). Positive real number. Can be a single number (all fleets) or a vector nfleets long (8).
#' @param Smax Maximum size limit (mm). Positive real number larger than Smin (if specified). Can be a single number (all fleets) or a vector nfleets long (8).
#' @param BL Bag limit (fish retained per trip). Positive integer. Can be a single number (all fleets) or a vector nfleets long (8). Note that predictive models are only available for RecS and HireS in regions FLK and NCFL.
#' @param PRM Post Release Mortality (discard mortality). Fraction. Can be a single number or a vector of length nareas (7).
#' @examples
#' mahiMP(1, Example_data)
#' @author T. Carruthers
#' @export
mahiMP = function(x, DataList,
                  Effort = NA,  TAC = NA,
                  Smin = NA, Smax = NA,
                  BL = NA, PRM = 0.85, ...){

  # Fleets (n=8): (1) USCom, (2) RecN, (3) RecS, (4) HireN, (5) HireS, (6) Intl, (7) Disc, (8) UnRep
  # Areas  (n=7): (1) NED, (2) NE, (3) NC, (4) SE, (5) SFL, (6) SAR, (7) CAR

  Finfo = EffDist #readRDS("C:/GitHub/DolphinMSE/MP_design/Finfo.rds")

  F_qfr = Finfo$F_qfr
  C_qfr = Finfo$C_qfr
  S_fa = Finfo$S_fa
  BGmu = Finfo$BGmu; BGcv = Finfo$BGcv

  qq = getq(DataList) # the next quarter

  mahiMPcheck(DataList, Effort, TAC, Smin, Smax, BL, F_qfr, S_fa) # check dimensions
  Rec = Blank_rec(DataList)                                       # make blank Rec object

  # no effort or TAC limits
  if(all(is.na(Effort[1]), is.na(TAC[1]))){
    cat("You did not specify any of these MP arguments: Effort, TAC.  Constant current F assumed \n")
    Rec = Rec_Effort(Rec, x, DataList, Effort = 1, BL, Smin, Smax, F_qfr, S_fa, BGmu, BGcv, qq)                        # Effort = 1 is current F essentially
  }

  # do Effort controls
  if(!is.na(Effort[1])){
    Rec = Rec_Effort(Rec, x, DataList, Effort, BL, Smin, Smax, F_qfr, S_fa, BGmu, BGcv, qq)        # Effort can be length 1, nf (8), nr (7), or c(nf,nr)
  }

  # do TAC controls
  if(!is.na(TAC[1])){
    Rec = Rec_TAC(Rec, x, DataList, TAC, BL, Smin, Smax, C_qfr, F_qfr, S_fa, BGmu, BGcv, qq)                    # TAC can be length 1, nf (8), nr (7), or c(nf,nr)
  }

  # Discard mortality rate
  Rec = Rec_Fdisc(Rec, PRM)

  # do Vessel Lims
  # only fleets 2-5 (RecN, RecS, HireN, HireS)
  Rec                                                                                  # return the hierarchical recommendation list, stock then fleet

}

class(mahiMP) = 'MMP'
