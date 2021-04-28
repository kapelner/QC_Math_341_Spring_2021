#problem 3
sum_x = 6
n = 20
xbar = sum_x / n
xbar

#problem 4
x = c(79,  91, 103,  87,  86,  91, 75)
mean(x)
n = length(x)
sigsq = 15

post_pred_sd = sqrt(sigsq + sigsq / n)
post_pred_sd

#problem 5
x = c(18,-3, 4)
n = length(x)
theta = 7.5
mean(x)
sigsqhatmle = sum((x - theta)^2) / n
sigsqhatmle
sigsqhatmle * n / 2
ssq = sum((x - mean(x))^2) / n
ssq * n / 2

#problem 6
alpha = 5.93
beta = 18.17
beta / (alpha - 1)
n0 = alpha * 2
sigsq0 = beta * 2 / n0
n0
sigsq0

(sigsqhatmle * n + sigsq0 * n0) / (n + n0 - 2)
(n0 - 2) / (n + n0 - 2)
