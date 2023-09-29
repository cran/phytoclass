## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(phytoclass)

## -----------------------------------------------------------------------------
Cluster.result <- Cluster(Sm, 14)

## ----fig.width=7--------------------------------------------------------------
# list of clusters
Cluster.result$cluster.list
# plot of clusters
plot(Cluster.result$cluster.plot)

## ----message=FALSE------------------------------------------------------------
set.seed("7683")
Results <- simulated_annealing(Sm, niter = 1)

## ----results------------------------------------------------------------------
Results$`condition number`
Results$RMSE
Results$MAE
Results$Error

Results$`F matrix`
Results$`Class abundances`

## ----figure-results, fig.width=7----------------------------------------------
Results$Figure

## ----message=FALSE------------------------------------------------------------
Clust1 <- Cluster(Sm, min_cluster_size = 14)$cluster.list[[1]]
# Remove the cluster column/label
Clust1$Clust <- NULL

set.seed("7683")
Results <- simulated_annealing(Clust1, niter = 1)

## ----results-clustering-------------------------------------------------------
Results$`condition number`
Results$RMSE
Results$MAE
Results$Error

Results$`F matrix`
Results$`Class abundances`

## ----figure-results-clustering, fig.width=7-----------------------------------
Results$Figure

## -----------------------------------------------------------------------------
#Create Fm (F matrix). Alternatively, a .csv file can be uploaded.
Fu <- data.frame(
  Per = c(0, 0, 0, 0, 1, 0, 0, 0),
  X19but = c(0, 0, 0, 0, 0, 1, 1, 0),
  Fuco = c(0, 0, 0, 1, 0, 1, 1, 0),
  Pra = c(1, 0, 0, 0, 0, 0, 0, 0),
  X19hex = c(0, 0, 0, 0, 0, 1, 0, 0),
  Allo = c(0, 0, 1, 0, 0, 0, 0, 0),
  Zea = c(1, 1, 0, 0, 0, 0, 0, 1),
  Chl.b = c(1, 1, 0, 0, 0, 0, 0, 0),
  Tchla = c(1, 1, 1, 1, 1, 1, 1, 1)
)

rownames(Fu) <- c(
  "Prasinophytes", "Chlorophytes", "Cryptophytes"
  , "Diatoms-B", "Dinoflagellates-A",
  "Haptophytes", "Pelagophytes", "Syn"
)

Min_max <- data.frame(
  Class = c(
    "Syn", "Chlorophytes", "Chlorophytes", "Prasinophytes", "Prasinophytes",
    "Prasinophytes", "Cryptophytes", "Diatoms-B", "Diatoms-B", "Pelagophytes",
    "Pelagophytes", "Pelagophytes", "Dinoflagellates-A", "Haptophytes",
    "Haptophytes", "Haptophytes", "Haptophytes", "Diatoms-B", "Cryptophytes",
    "Prasinophytes", "Chlorophytes", "Syn", "Dinoflagellates-A", "Pelagophytes"
  ),
  Pig_Abbrev = c(
    "Zea", "Zea", "Chl.b", "Pra", "Zea", "Chl.b", "Allo", "Chl.c3",
    "Fuco", "Chl.c3", "X19but", "Fuco", "Per", "X19but", "X19hex",
    "Fuco", "Tchla", "Tchla", "Tchla", "Tchla", "Tchla", "Tchla", "Tchla",
    "Tchla"
  ),
  min = as.numeric(c(
    0.0800, 0.0063, 0.1666, 0.0642, 0.0151, 0.4993, 0.2118, 0.0189,
    0.3315, 0.1471, 0.2457, 0.3092, 0.3421, 0.0819, 0.2107, 0.0090,
    1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000
  )),
  max = as.numeric(c(
    1.2123, 0.0722, 0.9254, 0.4369, 0.1396, 0.9072, 0.5479, 0.1840,
    0.9332, 0.2967, 1.0339, 1.2366, 0.8650, 0.2872, 1.3766, 0.4689,
    1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000
  ))
)

## ----message=FALSE------------------------------------------------------------
set.seed("7683")
Results <- simulated_annealing(
  S = Sm, 
  F = Fu,
  user_defined_min_max = Min_max,
  do_matrix_checks = TRUE,
  niter = 1,
  step = 0.01,
  weight.upper.bound = 30
)

## ----results-not-default------------------------------------------------------
Results$`condition number`
Results$RMSE
Results$MAE
Results$Error

Results$`F matrix`
Results$`Class abundances`

## ----figure-results-not-default, fig.width=7----------------------------------
Results$Figure

## -----------------------------------------------------------------------------
MC <- Matrix_checks(Sm, Fm)  
Snew <- MC$Snew
Fnew <- MC$Fnew  

## ----message=FALSE------------------------------------------------------------
MC <- Matrix_checks(Sm, Fm)  
Snew <- MC$Snew
Fnew <- MC$Fnew
SDRes <- Steepest_Desc(Fnew, Snew, num.loops = 10)

## -----------------------------------------------------------------------------
Bounded_weights(Sm, weight.upper.bound = 30)

## -----------------------------------------------------------------------------
MC <- Matrix_checks(Sm, Fm)  
Snew <- MC$Snew
Fnew <- MC$Fnew
cm <- Bounded_weights(Snew, weight.upper.bound = 30)

Results <- NNLS_MF(Fnew, Snew, cm)
Results$`F matrix`
Results$RMSE
Results$`C matrix`

