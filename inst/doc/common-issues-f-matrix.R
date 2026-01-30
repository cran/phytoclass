## -----------------------------------------------------------------------------
# # === remove numeric rownames introduced by read.csv
# if (all(grepl("^[0-9]+$", rownames(F_matrix)))) {
#   print("dropping unneeded numeric index")
#   # Set the first column as row names
#   rownames(F_matrix) <- F_matrix[[1]]
# 
#   # Remove the first column
#   F_matrix <- F_matrix[, -1]
# }

