#!/usr/bin/env Rscript

# Load libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stats)

# MODIFY PARAMETERS AS NEEDED
alpha <- 0.05
fnam <- "test_data.csv"
dir_path <- getwd()


#FUNCTION DECLARATIONS 
clear_screen <- function() {
  os <- Sys.info()["sysname"]
  
  if (os == "Windows") {
    system("cls")
  } else if (os == "Linux" || os == "Darwin") {  # Darwin for macOS
    cat("\014")
  } else {
    cat("Unsupported OS. Cannot clear screen.\n")
  }
}


fun <- function(x, t, g, c) {
  th <- exp(-x*t)
  sum(th * g) - c
}


# BEGIN CALCULATIONS
Z <- qnorm(1-alpha/2)

# Read data
df <- read_csv(fnam, col_types = cols())
t <- df$t
n <- df$n
h <- df$h

# Calculate variables
N <- n[1]
K <- sum(h)
g <- h/K
f <- c(n[-length(n)] - n[-1], 0)/N

# Create new data frame
my_data <- data.frame(t, n, h, g, f)
fnam2 <- paste0(gsub(".csv", "", fnam), "_added.csv")
write.csv(my_data, file = fnam2, row.names = FALSE)

# Calculate R0, Longevity, Generation time, and Population growth rate
R0 <- K/N
cat("R0 :", R0, "\n")

EL <- sum(t*f)
EL2 <- sum(t^2*f)
VL <- EL2 - EL^2
CI_L <- c(EL , VL, EL-Z*sqrt(VL/N), EL+Z*sqrt(VL/N))
cat("Longevity :", CI_L, "\n")

ET <- sum(t*g)
ET2 <- sum(t^2*g)
VT <- ET2 - ET^2
CI_mu <- c(ET, VT, ET-Z*sqrt(VT/K), ET+Z*sqrt(VT/K))
cat("Generation time :", CI_mu, "\n")

x_initial_guess <- 0.05  # initial guess for x
r_sol <- uniroot(function(x) fun(x, t, g, 1/R0), interval = c(0, 1), tol = 1e-6)$root

mu <- 1/R0
z <- exp(-(2*r_sol)*t)
s2 <- sum(z * g) - mu^2
s <- sqrt(s2)

r_L_sol <- uniroot(function(x) fun(x, t, g, mu+Z*s/sqrt(K)), interval = c(0, 1), tol = 1e-6)$root
r_U_sol <- uniroot(function(x) fun(x, t, g, mu-Z*s/sqrt(K)), interval = c(0, 1), tol = 1e-6)$root

CI_r <- c(r_sol, r_L_sol, r_U_sol)
cat("r :", CI_r, "\n")

CI_lam <- exp(CI_r)
cat("lambda:", CI_lam, "\n")

fun <- function(x, t, g, c) {
  th <- exp(-x*t)
  sum(th * g) - c
}

clear_screen()

# RESULTS
cat("Initial number of individuals N :", N, "\n")
cat("Offspring size K :", K, "\n")
cat("R0 :", R0, "\n")
cat("Longevity :", CI_L, "\n")
cat("Generation time :", CI_mu, "\n")
cat("r :", CI_r, "\n")
cat("lambda:", CI_lam, "\n")

cat("New data saved to: ", fnam2, "\n")


           