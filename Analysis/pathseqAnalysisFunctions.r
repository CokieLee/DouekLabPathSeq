library(vegan)

scaledPercentile <- function(countsCol) {
  n <- length(countsCol)
  rank <- rank(countsCol)
  percentile <- (rank - 1) / (n - 1)
  scaled_percentile <- percentile * (n / (n + 1))
  
  return(scaled_percentile)
}

basePlot <- function(title, xaxis, yaxis) {
  plot <-
    ggplot() +
    labs(title = title,
         x = xaxis,
         y = yaxis) +
    theme_bw()
  
  return(plot)
}