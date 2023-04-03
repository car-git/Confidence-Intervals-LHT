#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stats)

# Modify this as needed :---------------
alpha <- 0.05 # Change parameter as needed
fnam <- "test_data.csv" # Change file name
file <- fnam
# Change directory path here, live commented if file is in current directory:
dir_path <- getwd()   # use this if data is in current directory
#file <- file.path(dir_path, fnam)  # ....or change directory here

# Begin calculations :---------------

Z <- qnorm(1-alpha/2)
cat(Z, "\n")

fun <- function(x, t, g, c) {
  th <- exp(-x*t)
  sum(th * g) - c
}




df <- read_csv(file, col_types = cols())

t <- df$t
n <- df$n
h <- df$h

N <- n[1]
K <- sum(h)
g <- h/K

f <- c(n[-length(n)] - n[-1], 0)/N

my_data <- data.frame(t,n,h,g,f)
fnam2 <- paste0(gsub(".csv", "", fnam), "_added.csv")
write.csv(my_data, file = fnam2, row.names = FALSE)


# R0
R0 <- K/N
cat("R0 :", R0, "\n")

# Longevity 
EL <- sum(t*f)
EL2 <- sum(t^2*f)
VL <- EL2 - EL^2
CI_L <- c(EL , VL, EL-Z*sqrt(VL/N), EL+Z*sqrt(VL/N))
cat("Longevity :", CI_L, "\n")

# Generation time
ET <- sum(t*g)
ET2 <- sum(t^2*g)
VT <- ET2 - ET^2
CI_mu <- c(ET, VT, ET-Z*sqrt(VT/K), ET+Z*sqrt(VT/K))
cat("Generation time :", CI_mu, "\n")

# Population growth rate, r
x_initial_guess <- 0.05  # initial guess for x
r_sol <- uniroot(function(x) fun(x, t, g, 1/R0), interval = c(0, 1), tol = 1e-6)$root

mu <- 1/R0
z <- exp(-(2*r_sol)*t)
s2 <- sum(z * g) - mu^2
s <- sqrt(s2)

r_L_sol <- uniroot(function(x) fun(x, t, g, mu+Z*s/sqrt(K)), interval = c(0, 1), tol = 1e-6)$root
r_U_sol <- uniroot(function(x) fun(x, t, g, mu-Z*s/sqrt(K)), interval = c(0, 1), tol = 1e-6)$root

CI_r<- c(r_sol, r_L_sol, r_U_sol)
cat("r :", CI_r, "\n")

#CI_lam <- 1/thetas2
CI_lam <- exp(CI_r)
cat("lambda:", CI_lam, "\n")

cat("\014")   # Clear screen

# Results:
cat("Initial number of individuals N :", N)
cat("Offspring size K :", K)
cat("R0 :", R0)
cat("Longevity :",CI_L)
cat("Generation time :",CI_mu)
cat("r :", CI_r)
cat("lambda:", CI_lam)

cat("New data saved to:  ", fnam2, "\n")


