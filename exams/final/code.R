set.seed(1900)
n = 30

true_theta_0 = 0.1
true_theta_1 = 0.8
true_m = 7


x = rbinom(true_m, 1, true_theta_0)
x = c(x, rbinom(n - true_m, 1, true_theta_1))
sum_x = sum(x)

par(mfrow = c(1, 1))
plot(1 : n, x, type = "o")
#plot

## hyperparams
alpha = 1
beta = 1

#chains
S = 100000
theta_0s = array(NA, S)
theta_1s = array(NA, S)
ms = array(NA, S)
#start positions
theta_0s[1] = 0.5
theta_1s[1] = 0.5
ms[1] = 1

for (t in 2 : S){
  m = ms[t - 1]
  a_t = ifelse(m == 0, 0, sum(x[1 : m]))
  b_t = ifelse(m == n, 0, sum_x - a_t)
  
  theta_0 = rbeta(1, a_t + alpha, m - a_t + beta)
  theta_1 = rbeta(1, b_t + beta, n - m - b_t + beta)
  
  #now we need to calculate the m dist
  ln_m_dist = array(NA, n + 1)
  for (m in 0 : n){
    ln_m_dist[m + 1] = lchoose(m, a_t) + a_t * log(theta_0) + (m - a_t) * log(1 - theta_0) + 
      lchoose(n - m, b_t) + b_t * log(theta_1) + (n - m - b_t) * log(1 - theta_1)
    if (is.na(ln_m_dist[m + 1])){
      stop("boom")
    }
    
  }
  m_dist = exp(ln_m_dist) / sum(exp(ln_m_dist))
  
  theta_0s[t] = theta_0
  theta_1s[t] = theta_1
  ms[t] = sample(0 : n, 1, prob = m_dist)
}



###assess convergence
par(mfrow = c(3, 1))
S0 = 1000
plot(1 : S0, theta_0s[1 : S0])
# abline(h = mean(theta_0s[B : S0]), col = "blue")
# abline(h = true_theta_0, col = "red")
# abline(v = B, col = "grey")

plot(1 : S0, theta_1s[1 : S0])
# abline(h = mean(theta_1s[B : S0]), col = "blue")
# abline(h = true_theta_1, col = "red")
# abline(v = B, col = "grey")

plot(1 : S0, ms[1 : S0])
# abline(h = mean(ms[B : S0]), col = "blue")
# abline(h = sqrt(true_m), col = "red")
# abline(v = B, col = "grey")
#plot
BURN = 10

##assess autocorrelation

par(mfrow = c(3, 1))
xmax = 50
rho_max = 0.7
acf(theta_0s[BURN : S], xlim = c(1, xmax), ylim = c(-0.05, rho_max))
acf(theta_1s[BURN : S], xlim = c(1, xmax), ylim = c(-0.05, rho_max))
acf(ms[BURN : S], xlim = c(1, xmax), ylim = c(-0.05, rho_max))

THIN = 9
#plot

#burn and thin
theta_0s = theta_0s[(BURN + 1) : S]
theta_0s = theta_0s[seq(1, S - (BURN + 1), by = THIN)]
theta_1s = theta_1s[(BURN + 1) : S]
theta_1s = theta_1s[seq(1, S - (BURN + 1), by = THIN)]
ms = ms[(BURN + 1) : S]
ms = ms[seq(1, S - BURN, by = THIN)]


#look at posteriors with post-exp at 95% CI
par(mfrow = c(2, 1))


hist(theta_0s, br = 500, xlim = c(0, 1), xlab = "")
# abline(v = mean(theta_0s), col = "blue", lwd = 3)
# abline(v = quantile(theta_0s, 0.025), col = "grey", lwd = 3)
# abline(v = quantile(theta_0s, 0.975), col = "grey", lwd = 3)
# abline(v = true_theta_0, col = "red", lwd = 3)

hist(theta_1s, br = 500, xlim = c(0, 1), xlab = "")
# abline(v = mean(theta_1s), col = "blue", lwd = 3)
# abline(v = quantile(theta_1s, 0.025), col = "grey", lwd = 3)
# abline(v = quantile(theta_1s, 0.975), col = "grey", lwd = 3)
# abline(v = true_theta_1, col = "red", lwd = 3)

median(theta_1s)
mean(theta_1s)

par(mfrow = c(1, 1))
hist(ms, br = 500, xlim = c(0, n), xlab = "")
# abline(v = mean(ms), col = "blue", lwd = 3)
# abline(v = quantile(ms, 0.025), col = "grey", lwd = 3)
# abline(v = quantile(ms, 0.975), col = "grey", lwd = 3)
# abline(v = true_m, col = "red", lwd = 3)
#plot


par(mfrow = c(1, 1))
hist(theta_0s / theta_1s, br = 500, xlab = "")












####################
rm(list = ls())

pacman::p_load(MCMCpack, ggplot2)

set.seed(1)
n = 45
x = rnorm(n, .11, .06)
xbar = mean(x)
s_sq = var(x)

S = 10000
sigsq_samps = rinvgamma(S, n / 2, (n-1) * s_sq / 2)
theta_samps = rnorm(S, xbar, sqrt(sigsq_samps / n))

sum(sigsq_samps < 0.045^2) / S
ggplot(data.frame(theta = theta_samps, sigsq = sigsq_samps)[1:1000,]) +
  aes(x = theta, y = sigsq) + 
  geom_point()

ggplot(data.frame(theta = theta_samps, sigsq = sigsq_samps)) +
  aes(x = theta, y = sigsq) + 
  geom_point() + 
  ylim(0.002, 0.003) + 
  geom_smooth(se = FALSE, lwd =2)
  # geom_abline(slope = coef(mod)[2], intercept = coef(mod)[1],col = "blue", lwd = 2)


# set.seed(1)
# n = 16
# x = rnorm(n, .11, .06)
# xbar = mean(x)
# s = sd(x)

