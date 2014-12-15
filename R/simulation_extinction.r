# exponential survival
# S(t) = exp^{-lambda t}
# h(t) = lambda
exp.surv <- function(time, lambda) {
  exp(-(lambda) * time)
}

exp.haz <- function(time, lambda) {
  lambda
}


# Weibull survival
# S(t) = exp^{-lambda t^p}
# h(t) = lambda p t^{p - 1}
wei.surv <- function(time, lambda, p) {
  tt <- time**p
  exp(-(lambda) * tt)
}

wei.haz <- function(time, lambda, p) {
  tt <- time**(p - 1)
  lambda * p * tt
}


rng <- data.frame(x <- c(0, 12))
xx <- seq(0, 5, by = 0.01)
wz <- wei.haz(xx, 1, 0.9)
plot(xx, wz)
