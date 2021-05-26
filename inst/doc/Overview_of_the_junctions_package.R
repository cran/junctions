## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.width = 7, fig.height = 7, echo = TRUE)
library(junctions)
require(tidyr)
require(dplyr)
require(ggplot2)
require(magrittr)

## ----phased unphased simulation-----------------------------------------------
simulated_pop <- sim_phased_unphased(pop_size = 1000,
                                     freq_ancestor_1 = 0.5,
                                     total_runtime = 100,
                                     size_in_morgan = 1,
                                     markers = sort(runif(n = 1000, 0, 1)),
                                     time_points = seq(0, 100, by = 10))

## ----sim fin chrom------------------------------------------------------------
junctions_fin_chrom <- sim_fin_chrom(pop_size = 1000,
                                     freq_ancestor_1 = 0.5,
                                     total_runtime = 200,
                                     morgan = 1,
                                     R = 100) # the number of crossover sites

junctions_inf_chrom <- sim_inf_chrom(pop_size = 1000,
                                     freq_ancestor_1 = 0.5,
                                     total_runtime = 200,
                                     morgan = 1,
                                     markers = 100)

plot(junctions_inf_chrom$avgJunctions,
     type = "l", lwd = 2,
     col = "darkgreen",
     xlab = "Time since admixture",
     ylab = "Number of junctions")
lines(junctions_fin_chrom$avgJunctions, col = "blue", lwd = 2)
lines(junctions_inf_chrom$detectedJunctions, col = "lightgreen", lwd = 2)
legend("topleft", legend = c("Finite chromosome", "Infinite chromosome",
                             "Infinite chromosome detected"),
       col = c("blue", "darkgreen", "lightgreen"), lty = 1, lwd = 2)

## ----backcrossing-------------------------------------------------------------
backcross_result <- sim_backcrossing(population_size = 1000,
                                     freq_ancestor_1 = 0.5,
                                     total_runtime = 20,
                                     size_in_morgan = 1,
                                     number_of_markers = 100)

plot(backcross_result$average_junctions,
     type = "l", lwd = 2,
     col = "darkgreen",
     xlab = "Time since admixture",
     ylab = "Number of junctions")

## ----estimate_time------------------------------------------------------------
estimate_time(J = 10, N = 1000, R = 1000, H_0 = 0.5, C = 1)

## ----estimate time markers----------------------------------------------------
estimate_time_one_chrom(J = 10, N = 1000, H_0 = 0.5,
                        marker_distribution = sort(runif(n = 1000, 0, 1)))

## ----estimate time haploid----------------------------------------------------
ancestry_data <- subset(simulated_pop, simulated_pop$time == 100)
ancestry_matrix <- dplyr::select(ancestry_data, c(individual,
                                                  location, anc_chrom_1))
estimate_time_haploid(ancestry_matrix = ancestry_matrix,
                      N = 1000,
                      freq_ancestor_1 = 0.5)

## ----estimate time diploid----------------------------------------------------
ancestry_matrix <- dplyr::select(ancestry_data, c(individual,
                                                  location,
                                                  anc_chrom_1, anc_chrom_2))
ancestry_matrix <- cbind(rep(1, length(ancestry_matrix$individual)),
                         ancestry_matrix)

t_phased <- estimate_time_diploid(ancestry_information = ancestry_matrix,
                                  analysis_type = "individuals",
                                  phased = TRUE,
                                  pop_size = 1000,
                                  freq_ancestor_1 = 0.5)

t_unphased <- estimate_time_diploid(ancestry_information = ancestry_matrix,
                                    analysis_type = "individuals",
                                    phased = FALSE,
                                    pop_size = 1000,
                                    freq_ancestor_1 = 0.5)

t_phased
t_unphased

## ----likelihood---------------------------------------------------------------
ancestry_matrix <- dplyr::select(ancestry_data, c(individual, location,
                                                  anc_chrom_1, anc_chrom_2))

time_points <- 80:120
ll_phased <- log_likelihood_diploid(ancestry_matrix,
                                    pop_size = 1000,
                                    freq_ancestor_1 = 0.5,
                                    t = time_points,
                                    phased = TRUE)
ll_unphased <- log_likelihood_diploid(ancestry_matrix,
                                      pop_size = 1000,
                                      freq_ancestor_1 = 0.5,
                                      t = time_points,
                                      phased = FALSE)

to_plot <- tibble::tibble(time_points, ll_phased, ll_unphased)
to_plot %>%
  tidyr::gather(key = "phasing", value = "loglikelihood", -time_points) %>%
  ggplot2::ggplot(
      ggplot2::aes(x = time_points, y = loglikelihood, col = phasing)) +
  ggplot2::geom_line()

## ----likelihood haploid-------------------------------------------------------
ancestry_matrix <- dplyr::select(ancestry_data, c(individual, location,
                                                  anc_chrom_1))

ll_haploid <- log_likelihood_haploid(ancestry_matrix,
                                     N = 1000,
                                     freq_ancestor_1 = 0.5,
                                     t = time_points)
plot(ll_haploid ~ time_points, type = "l",
     xlab = "Time since admixture",
     ylab = "Loglikelihood")

## ----expected number of junctions---------------------------------------------
number_of_junctions(N = 100, R = 100, H_0 = 0.5, C = 1, t = 1000)

## ----expected number of junctions markers-------------------------------------
number_of_junctions_markers(N = 100, H_0 = 0.5, t = 1000,
                            marker_distribution = sort(runif(100, 0, 1)))

## ----expected number of junctions di------------------------------------------
number_of_junctions_di(N = 100, H_0 = 0.5, t = 1000, di = 1e-5)

## ----expected number of junctions back crossing-------------------------------
number_of_junctions_backcross(H_0 = 0.5, C = 1, t = 10)

## ----K------------------------------------------------------------------------
calc_k(N = 1000, R = 1000, H_0 = 0.5, C = 1)

## ----MAT----------------------------------------------------------------------
calculate_mat(N = 1000, R = 1000, H_0 = 0.5, C = 1)

## ----error--------------------------------------------------------------------
time_error(t = 30, N = 1000, R = 1000, H_0 = 0.5, C = 1)

