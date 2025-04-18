## ----setup, include=FALSE-----------------------------------------------------
library(junctions)
knitr::opts_chunk$set(echo = TRUE)

## ----eval = FALSE-------------------------------------------------------------
# sim_inf_chrom(pop_size,
#               initial_heterozygosity,
#               total_runtime,
#               morgan,
#               markers,
#               seed)

## ----single_sim_inf, fig.width = 7, fig.height = 7----------------------------
pop_size <- 100 # population size
h_0 <- 0.5 # initial heterozygosity
maximum_time <- 1000 # run time
c <- 1 # number of recombinations per meiosis

# we first ignore the effect of imposing randomly distributed markers
number_of_markers <- -1

v <- sim_inf_chrom(pop_size = pop_size,
                   freq_ancestor_1 = h_0,
                   total_runtime = maximum_time,
                   morgan = c,
                   markers = number_of_markers,
                   seed = 42)

plot(v$avgJunctions,
     type = "l",
     xlab = "Generations",
     ylab = "Number of Junctions",
     main = "Example Infinite Chromosome")

## ----repl_sim_inf-------------------------------------------------------------
number_replicates <- 10
v <- c()
for (r in 1:number_replicates) {
  v2 <- sim_inf_chrom(pop_size = pop_size,
                      freq_ancestor_1 = h_0,
                      total_runtime = maximum_time,
                      morgan = c,
                      markers = number_of_markers,
                      seed = r)
  v <- rbind(v, as.numeric(v2$avgJunctions))
}
v <- colMeans(v) #mean across replicates

## ----plot_sim_inf, fig.width = 7, fig.height = 7------------------------------
clarity <- seq(1,
               maximum_time,
               length.out = 50) # we plot not all points, for clarity
plot(v[clarity] ~ clarity, lwd = 2,
     xlab = "Generations",
     ylab = "Number of Junctions",
     main = "Average behaviour Infinite chromosome", pch = 16)

t <- 0:maximum_time
predicted <- number_of_junctions(N = pop_size, H_0 = h_0, C = c, t = t)
lines(predicted ~ t, col = "blue")

## ----single_sim_random, fig.width = 7, fig.height = 7-------------------------
pop_size <- 100 # population size
h_0 <- 0.5 # initial heterozygosity
maximum_time <- 1000 # run time
c <- 1 # number of recombinations per meiosis
number_of_markers <- 1000 # 1000 markers

#single example run
v <- sim_inf_chrom(pop_size = pop_size,
                   freq_ancestor_1 = h_0,
                   total_runtime = maximum_time,
                   morgan = c,
                   markers = number_of_markers,
                   seed = 42)
plot(v$avgJunctions,
     type = "l",
     xlab = "Generations",
     ylab = "Number of Junctions",
     main = "Example Infinite Chromosome")

lines(v$detectedJunctions, col = "blue")
legend("bottomright",
       c("Real number", "Number detected"),
       lty = 1,
       col = c("black", "blue"))

## ----repl_sim_random----------------------------------------------------------
mean_junctions <- c()
detected_junctions <- c()
for (r in 1:number_replicates) {
  v2 <- sim_inf_chrom(pop_size = pop_size,
                      freq_ancestor_1 =  h_0,
                      total_runtime = maximum_time,
                      morgan = c,
                      markers = number_of_markers,
                      seed = r + 42)
  mean_junctions <- rbind(mean_junctions,
                          as.numeric(v2$avgJunctions))
  detected_junctions <- rbind(detected_junctions,
                              as.numeric(v2$detectedJunctions))
}
mean_junctions <- colMeans(mean_junctions)
detected_junctions <- colMeans(detected_junctions)

## ----plot_sim_random, fig.width = 7, fig.height = 7---------------------------
#we plot not all points, for clarity
clarity <- seq(1, maximum_time, length.out = 50)
plot(mean_junctions[clarity] ~ clarity,
     xlab = "Generations",
     ylab = "Number of Junctions",
     main = "Average behaviour Infinite chromosome",
     pch = 16)

points(detected_junctions[clarity] ~ clarity,
       pch = 17,
       col = "blue")

t <- 0:maximum_time
predicted <- number_of_junctions(N = pop_size,
                                 H_0 = h_0,
                                 C = c,
                                 t = t)
lines(predicted ~ t, lwd = 2)

# now substitute K with the observed maximum detected.
k <- tail(detected_junctions, 1)
pred <- k - k * (1 - h_0 * c / k) ^ t
lines(pred ~ t, col = "blue", lwd = 2)

legend("bottomright",
       c("Real number", "Detected",
         "Predicted Real", "Predicted Detected"),
       pch = c(16, 17, NA, NA),
       lty = c(NA, NA, 1, 1),
       col = c("black", "blue",
               "black", "blue"),
       lwd = 2)

## ----eval = FALSE-------------------------------------------------------------
# sim_fin_chrom(pop_size,
#               initial_heterozygosity,
#               total_runtime,
#               morgan,
#               seed,
#               R)

## ----single_sim_fin, fig.width = 7, fig.height = 7----------------------------
r <- 100 # chromosome size
n <- 100 # population size
freq_ancestor_1 <- 0.5 # frequency of ancestor 1 at t = 0
c <- 1 # number of recombinations per meiosis
maximum_time <- 1000

#single example run
v <- sim_fin_chrom(pop_size = n,
                   freq_ancestor_1 = freq_ancestor_1,
                   total_runtime = maximum_time,
                   morgan = c,
                   seed = 42,
                   R = r)
plot(v$avgJunctions, type = "l",
     xlab = "Generations",
     ylab = "Number of Junctions",
     main = "Example Finite Chromosome")

## ----repl_sim_fin-------------------------------------------------------------
v <- c()
for (repl in 1:number_replicates) {
  v2 <- sim_fin_chrom(pop_size = n,
                      freq_ancestor_1 = h_0,
                      total_runtime = maximum_time,
                      morgan = c,
                      seed = repl,
                      R = r)
  v <- rbind(v, as.numeric(v2$avgJunctions))
}
v <- colMeans(v)

## ----plot_sim_fin, fig.width = 7, fig.height = 7------------------------------
clarity <- seq(1, 1000, length.out = 50) #we plot not all points, for clarity
plot(v[clarity] ~ clarity, lwd = 2,
     xlab = "Generations",
     ylab = "Number of Junctions",
     main = "Average behaviour Finite Chromosome",
     pch  = 16)

t <- 0:maximum_time
predicted <- number_of_junctions(N = n, R = r,
                                 H_0 = h_0, C = c,
                                 t)
lines(predicted ~ t, col = "blue")
legend("bottomright", c("Simulated", "Predicted"),
       pch = c(16, NA),
       lty = c(NA, 1),
       col = c("black", "blue"))

## ----equations----------------------------------------------------------------
num_j <- number_of_junctions(N = 100,
                             R = 1000,
                             H_0 = 0.5,
                             C = 1,
                             t = 200)
num_j <- tail(num_j, 1)
num_j

## -----------------------------------------------------------------------------
time_estim <- estimate_time(J = num_j,
                            N = 100,
                            R = 1000,
                            H_0 = 0.5,
                            C = 1)
time_estim

## -----------------------------------------------------------------------------
time_error(t = time_estim,
           N = 100,
           R = 1000,
           H_0 = 0.5,
           C = 1,
           relative = TRUE)

## -----------------------------------------------------------------------------
time_error(t = time_estim,
           N = 100,
           R = 1000,
           H_0 = 0.5,
           C = 1,
           relative = FALSE)

## -----------------------------------------------------------------------------
calculate_mat(N = 100,
              R = 1000,
              H_0 = 0.5,
              C = 1)

## ----fig.width = 7, fig.height = 7--------------------------------------------
maximum_time <- 1000
t <- 0:maximum_time
num_j <- number_of_junctions(N = 100,
                             R = 1000,
                             H_0 = 0.5,
                             C = 1,
                             t = t)
par(mar = c(4, 5, 2, 3))
plot(num_j ~ t,
     type = "l",
     xlab = "Generations",
     ylab = "Number of Junctions",
     xlim = c(0, maximum_time))
#vertical line that indicates the upper limit
abline(v = calculate_mat(N = 100,
                         R = 1000,
                         H_0 = 0.5,
                         C = 1),
       lty = 2)

par(new = TRUE)
v <- time_error(t = 0:(maximum_time - 1), #to avoid an error at t = maxT
                N = 100,
                R = 1000,
                H_0 = 0.5,
                C = 1,
                relative = TRUE)
plot(v,
     col = "red", type = "l",
     xlim = c(0, maximum_time),
     ylim = c(0, 1),
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "")
axis(4)
mtext("Relative error", side = 4, line = 2)
legend("bottomright",
       c("Number of junctions", "Relative Error", "t_MAX"),
       lty = c(1, 1, 2),
       col = c("black", "red", "black"))

