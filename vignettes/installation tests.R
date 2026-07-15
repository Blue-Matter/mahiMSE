


library(remotes)
options(timeout=500) # need to extend timeout or large package install won't work

install_github('blue-matter/MSEtool', ref='dev')               # dev branch of MSEtool
install_github("blue-matter/mahiMSE", build_vignettes = T)     # mahiMP, plotting etc
install_github("blue-matter/mahiRefSet")  # 32 OMs ~ 5 mins depending on connection
install_github("blue-matter/mahiRobSet")  # 13 OMs ~ 2 mins depending on connection

library(mahiMSE)
mahiMP(Example_Data)
Hist = Simulate(smallOM)
Proj = Project(Hist,'mahiMP')



vignette('mahiMSE')




Hist = Simulate(smallOM)

Proj = Project(Hist,'mahiMP')

mahiMP(Example_Data)
