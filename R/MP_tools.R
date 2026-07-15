



#' An MP that stores a data object while a projection is running
#'
#' @param Data An MSEtool object of class 'Data'
#' @param repyr The projected time step in which the data is written to disk. Positive integer.
#' @param outdir A character string representing the file directory to save Example_Data to.
#' @param filenam A character string representing the file name
#' @param Eff Real number, the constant effort for projection
#' @author T. Carruthers
#' @export
fakeMP = function(Data, repyr = 2023, outdir = "C:/GitHub/DolphinMSE/MP_design", filenam = "Example_Data", Eff=0.5){

  ad=Advice(Effort=Eff)
  if(max(Data@Years)==repyr){
    Example_Data = Data
    outfile =  paste0(outdir,"/",filenam,".rds")
    saveRDS(Example_Data, outfile)
    cat(paste0("Example data file with data in year ",repyr, " outputted to: ",outfile, "\n"))
    stop()
  }
  ad

}
class(fakeMP) = "mp"


#' MP helper function - what season is it?
#'
#' @param Data An MSEtool object of class 'Data'
#' @author T. Carruthers
#' @export
get_season = function(Data){
  rep(1:4,1000)[max(Data@Years)-Data@YearLH+1] # This is the season following the data of this season - so one time step ahead and the correct index for TAC and Effort distribution
}

#' A basic effort MP for testing
#'
#' @param Data An MSEtool object of class 'Data'
#' @param Eff Real number, the constant effort for projection
#' @author T. Carruthers
#' @export
testMP05 = function(Data, Eff=0.5){
 Advice(Effort=Eff)
}
class(testMP05) = "mp"

#' A basic effort MP for testing
#'
#' @param Data An MSEtool object of class 'Data'
#' @param Eff Real number, the constant effort for projection
#' @author T. Carruthers
#' @export
testMP01 = function(Data, Eff=0.1){
  Advice(Effort=Eff)
}
class(testMP01) = "mp"
