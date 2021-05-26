## ----setup, include=FALSE-----------------------------------------------------
library(junctions)
library(Rcpp)
knitr::opts_chunk$set(fig.width = 7, echo = TRUE)

## ---- sim_data----------------------------------------------------------------
simulated_data <- sim_phased_unphased(pop_size = 1000,
                                      freq_ancestor_1 = 0.5,
                                      total_runtime = 100,
                                      size_in_morgan = 1,
                                      time_points = 100,
                                      markers = 1000)

simulated_data

## ---- infer unphased admixture time-------------------------------------------
focal_data <- subset(simulated_data, simulated_data$time == 100)
admixture_time <- estimate_time_diploid(ancestry_information =
                                          cbind(1,
                                                1,
                                                focal_data$location,
                                                focal_data$anc_chrom_1,
                                                focal_data$anc_chrom_2),
                                         pop_size = 1000,
                                         freq_ancestor_1 = 0.5,
                                         phased = FALSE)
admixture_time

## ---- infer phased admixture time---------------------------------------------
morgan_locations <- focal_data$location
phased_data <- cbind(focal_data$anc_chrom_1, focal_data$anc_chrom_2)
admixture_time_phased <- estimate_time_diploid(ancestry_information =
                                                 cbind(1,
                                                       1,
                                                       morgan_locations,
                                                       phased_data),
                                              phased = TRUE,
                                              pop_size = 1000,
                                              freq_ancestor_1 = 0.5)

admixture_time_phased

## ---- infer time pop size-----------------------------------------------------
found <- c()
for (N in c(100, 1000, 10000, 100000, 1e6, 1e7)) {
  admixture_time_phased <- estimate_time_diploid(ancestry_information =
                                                 cbind(1,
                                                       1,
                                                       morgan_locations,
                                                       phased_data),
                                              phased = TRUE,
                                              pop_size = N,
                                              freq_ancestor_1 = 0.5)
  found <- rbind(found, c(N, admixture_time_phased$time[[1]],
                          admixture_time_phased$loglikelihood[[1]]))
}
found
plot(found[, 2] ~ found[, 1], log = "x",
     xlab = "Population Size",
     ylab = "Admixture Time")
plot((-1 * found[, 3]) ~ found[, 1], log = "x",
     xlab = "Population Size",
     ylab = "Log Likelihood")

## ---- likelihood--------------------------------------------------------------
found <- c()
for (N in 10 ^ (seq(1, 6, length.out = 100))) {
  ll <- junctions::log_likelihood_diploid(local_anc_matrix =
                                            cbind(1,
                                                  morgan_locations,
                                                  phased_data),
                                        phased = TRUE,
                                        pop_size = N,
                                        freq_ancestor_1 = 0.5,
                                        t = 100)
  found <- rbind(found, c(N, ll))
}
plot(found, xlab = "Population Size", ylab = "Log Likelihood", log = "x")

## ---- phasing error-----------------------------------------------------------
simulated_data <- sim_phased_unphased(pop_size = 1000,
                                        freq_ancestor_1 = 0.5,
                                        total_runtime = 100,
                                        size_in_morgan = 1,
                                        time_points = 100,
                                        markers = 1000,
                                        error_rate = 0.01)

simulated_data$true_data
simulated_data$phased_data

## ---- compare-----------------------------------------------------------------
focal_true_data <- subset(simulated_data$true_data,
                          simulated_data$true_data$individual == 0)

true_data <- cbind(focal_true_data$anc_chrom_1,
                   focal_true_data$anc_chrom_2)

true_loc  <- focal_true_data$location

admixture_time_true <- estimate_time_diploid(ancestry_information  =
                                               cbind(1,
                                                     1,
                                                     true_loc,
                                                     true_data),
                                            phased = TRUE,
                                            pop_size = 1000,
                                            freq_ancestor_1 = 0.5)

focal_phased_data <- subset(simulated_data$phased_data,
                            simulated_data$phased_data$individual == 0)

phased_data <- cbind(focal_phased_data$anc_chrom_1,
                     focal_phased_data$anc_chrom_2)
phased_loc  <- focal_phased_data$location


admixture_time_error <- estimate_time_diploid(ancestry_information =
                                               cbind(1,
                                                     1,
                                                     phased_loc,
                                                     phased_data),
                                            phased = TRUE,
                                            pop_size = 1000,
                                            freq_ancestor_1 = 0.5)
admixture_time_true
admixture_time_error

