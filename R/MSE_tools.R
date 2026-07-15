

# Export necessary code and objects
mahi_setup = function(MPs, cpus = 8){
  if(sfIsRunning())sfStop()
  library(snowfall)
  setup(cpus = cpus)
  sfLibrary(mahiRefSet)
  sfLibrary(mahiMSE)
  sfExport(list = MPs)
}

# Parallel projection
parallel_Project = function(X, MPs, largedir, nSim){
  OM = get(paste0("OM_",X))
  OM@Seed = X
  Hist = Simulate(OM, nSim=nSim)
  MSE = Project(Hist, MPs)
  saveRDS(MSE, paste0(largedir,"/MSE_",X,".rds"))
}
