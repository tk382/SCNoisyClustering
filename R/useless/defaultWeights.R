defaultWeights <- function (n)
{
  a <- seq(length = n - 1)/n
  b <- a * (1 - a)
  1/sqrt(b * n)
}