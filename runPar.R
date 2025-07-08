# Note: the following script assumes that the models has already been compiled 
# in a folder named "build" within the model folder.

# Load libraries
library(foreach)
library(doParallel)

# Set number of cores and start parallel backend
ncores <- detectCores()
cl <- makeCluster(ncores, outfile = "")
registerDoParallel(cl)

# ===== USER SETTINGS =====

# Modify below to specify the absolute path to the model folder
model_folder_path <- "[PATH]"

# Modify below to specify a different input file (assumed inside an input_file
# folder within the model folder)
inputfile <- "dsk_sfc_inputs_baseline.json"
input_file_path <- file.path(model_folder_path, "input_files", inputfile)

# Number of Monte Carlo simulations
MC <- 10

# Name for the run
runname <- "partest"

# ===== PARALLEL EXECUTION =====

foreach(j = 1:MC, .errorhandling = 'remove', .multicombine = TRUE) %dopar% {
  
  # Seed for each run
  seed <- paste("-s", j)
  
  # Path to executable
  exec <- file.path(model_folder_path, "build", "dsk_SFC_G")
  
  # Build system command
  command <- paste(exec, input_file_path, "-r", runname, seed)
  
  # Run the model
  system(command)
}

# Stop parallel backend
stopCluster(cl)