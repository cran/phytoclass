% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cluster.R
\name{Cluster}
\alias{Cluster}
\title{Cluster things}
\usage{
Cluster(Data, min_cluster_size)
}
\arguments{
\item{Data}{S (sample) matrix}

\item{min_cluster_size}{the minimum size required for a cluster}
}
\value{
A named list of length two. The first element "cluster.list"
is a list of clusters, and the second element "cluster.plot" the
cluster analysis object (dendogram) that can be plotted.
}
\description{
Cluster things
}
\examples{
Cluster.result <- Cluster(Sm, 14)
Cluster.result$cluster.list
plot(Cluster.result$cluster.plot)
}
