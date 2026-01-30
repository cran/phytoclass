## ----setup--------------------------------------------------------------------
library(phytoclass)

library(here)

## -----------------------------------------------------------------------------
# Load the sample matrix from a file
# If your file does not contain sample names (e.g., Sample_1, Sample_2), you may omit `row.names = 1`
# In that case, do: S_matrix <- read.csv(here("vignettes/custom-example-S.csv"))
S_matrix <- read.csv(here("vignettes/custom-example-S.csv"), row.names = 1)  

## -----------------------------------------------------------------------------
# pigment-taxa occurrence matrix (specific to an oceanographic region)
F_matrix <- read.csv(here("vignettes/custom-example-F.csv"))

## -----------------------------------------------------------------------------
# === remove numeric rownames introduced by read.csv
if (all(grepl("^[0-9]+$", rownames(F_matrix)))) {
  print("dropping unneeded numeric index")
  # Set the first column as row names
  rownames(F_matrix) <- F_matrix[[1]]
  
  # Remove the first column
  F_matrix <- F_matrix[, -1] 
}

## -----------------------------------------------------------------------------
min_max_matrix <- read.csv(here("vignettes/custom-example-min_max-southern_ocean.csv"))

## -----------------------------------------------------------------------------
# === remove numeric rownames introduced by read.csv
if (all(grepl("^[0-9]+$", rownames(min_max_matrix)))) {
  print("dropping unneeded numeric index")
  # Set the first column as row names
  rownames(min_max_matrix) <- min_max_matrix[[1]]
  
  # Remove the first column
  min_max_matrix <- min_max_matrix[, -1] 
}

## -----------------------------------------------------------------------------
phytoclass::simulated_annealing(
  S_matrix,
  Fmat = F_matrix,
  user_defined_min_max = min_max_matrix,
  verbose = FALSE
)

