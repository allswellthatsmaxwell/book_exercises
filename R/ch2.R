library(ggplot2)
library(magrittr)

inverse <- function(f, lower, upper){
  function(y){
    uniroot(function(x){f(x) - y}, lower = lower, upper = upper, tol=1e-3)[1]
  }
}

# ch2

# problem 1
# Posterior inference: suppose you have a Beta(4, 4) prior distribution on the probability θ
# that a coin will yield a ‘head’ when spun in a specified manner. The coin is independently
# spun ten times, and ‘heads’ appear fewer than 3 times. You are not told how many heads
# were seen, only that the number is less than 3. Calculate your exact posterior density
# (up to a proportionality constant) for θ and sketch it.
#
# Answer: Average the beta distributions for every possible world. The different worlds...
# might have different weights on them? I think they do? Since there are more ways to 
# flip three heads than ways to flip one head - but I don't want to think about that
# so I'm gonna skip it. 
xs <- seq(0, 1, 0.01)
priors <- list(a = 4, b = 4)
possible_heads <- 0:3
flips <- 10
possible_distributions <- lapply(
  possible_heads, 
  function(nheads) dbeta(xs, priors$a + nheads, 
                         priors$b + flips - nheads))
average_distribution <- (
  Reduce(function(d1, d2) d1 + d2, possible_distributions) / 
  length(possible_distributions))

plot(xs, possible_distributions[[1]], type = 'l', col = 'gray')
for (d in possible_distributions) {
  lines(xs, d, col = 'gray')
}
lines(xs, average_distribution, col = 'red')


# problem 2
# Predictive distributions: consider two coins, C1 and C2, with the following characteristics:
# Pr(heads|C1)=0.6 and Pr(heads|C2)=0.4. Choose one of the coins at random and
# imagine spinning it repeatedly. Given that the first two spins from the chosen coin are
# tails, what is the expectation of the number of additional spins until a head shows up?
#
# Answer: Two parts to this one. First: find the posterior after the two Tails spins.
# Second: compute the mean of the .. predictive distribution of that posterior, I think?
# P(C1|tails) = (P(tails|C1) * P(C1) / a) = 0.4 * 0.5 / a, where
#   a =([P(tails|C1) * P(C1) + P(tails|C2) * P(C2)]) = 0.4 * 0.5 + 0.6 * 0.5 = 
#   4/10 * 1/2 + 6/10 * 1/2 = 4/20 + 6/20 = 10/20 = 1/2, so we have 
#   the full expression (4/10 * 5/10) / 1/2 = (1/5) / (1/2) = 2/5.
# P(C1|T,T) = P(tails|C1) * P(C1|tails) / b = 4/10 * 2/5 / b, where b = 
# P(tails|C1) * P(C1|tails) + P(tails|C2) * P(C2|tails) = 
# 4/10 * 2/5 + 6/10 * 3/5 = 8/50 + 18/50 = 26/50 = b, so the full expression
# is 4/10 * 2/5 / 26/50 = 8/50 * 50/26 = 8/26 = 4/13!
# So now... we think it's a 4/13 chance that it's P(heads|all the above) = 0.6,
# and a 9/13 chance it's P(heads|all the above = 0.4). We can just weight these worlds,
# so our probability of heads on any flip is 4/13 * 6/10 + 9/13 * 4/10 = (24 + 36)/130 = 
# 60/130.
# So we want... the number of flips so ... oh!! use exponential distribution maybe! 
# the rate is 60/130!
# Or... this is where hierarchical modeling comes in. 
# We want to simulate a draw from the 4/13 vs. 9/13, then the draw from the proper exponential, 
# many times. Maybe?
outer_ndraws <- 100000
outer_draws <- runif(outer_ndraws) < 4/13 # TRUE means "picked C1"
coin_probas <- sapply(outer_draws, function(draw) if (draw) 0.6 else 0.4)
mean(sapply(coin_probas, function(r) rexp(1, r)))

# problem 3
times <- 1000
xs <- seq(0, 0.3, 0.001)
mu <- 1/6
variance <- sqrt(times * 1/6 * 5/6) / times
plot(xs, dnorm(xs, mu, variance), type='l')
x_locations <- c(0.05, 0.25, 0.50, 0.75, 0.95)
print(qnorm(x_locations, mu, variance))

# problem 4a
xs <- seq(0, 1, 0.001)
means <- c(1/12, 1/6, 1/4)
weights <- c(0.25, 0.5, 0.25)
distributions <- Map(function(mu, weight) weight * dnorm(xs, mu, variance) / times,
                     means, weights)
mixed <- Reduce(`+`, distributions)
plot(xs, mixed, type = 'l')

cumul_curve <- cumsum(mixed)
x_locations <- c(0.1, 0.2, 0.25, 0.3, 0.4)
quantiles <- quantile(cumul_curve, x_locations)
plot(xs, cumul_curve, type = 'l', xlim = c(0, 0.4))
for (i in seq_along(x_locations)) {
  x_loc <- x_locations[i]
  q <- quantiles[i]
  abline(v = x_loc, h = q, col = 'gray')
}

# problem 8c
compute_posterior <- function(n) {
  known_sd <- sigma <- 20
  prior_mean <- mu0 <- 180
  prior_sd <- tau0 <- 40
  data_mean <- ybar <- 150
  posterion_mean_denom <- (1 / tau0^2) + (n / sigma^2)
  posterion_mean_numer <- (mu0 / tau0^2) + ((ybar * n) / sigma^2)
  posterior_mean <- mu_n <- posterion_mean_numer / posterion_mean_denom
  
  posterior_sd_denom <- (1 / tau0^2) + (n / sigma^2)
  posterior_sd <- tau_n <- 1 / posterior_sd_denom
  list(mu = posterior_mean, sigma = posterior_sd)
}

get_posterior_interval <- function(n) {
  posterior <- compute_posterior(n)
  posterior %$% qnorm(c(0.05, 0.95), mu, sigma)
}

n <- 100
results <- list(posterior = compute_posterior(n), interval = get_posterior_interval(n))
xs <- seq(0, 350, 1)
distr <- results %$% posterior %$% dnorm(xs, mu, sigma)
plot(xs, distr, type = 'l')
abline(v=results$interval[[1]])
abline(v=results$interval[[2]])
results

# problem 9
alpha <- 1
beta <- 2/3
xs <- seq(0, 1, 0.0001)
dens <- dbeta(xs, alpha, beta)
points <- rbeta(1000, alpha, beta)
mean(points)
plot(xs, dens, type = 'l')

posterior <- list(a = 651, b = 350.67)
dens <- posterior %$% dbeta(xs, a, b)
rvec <- posterior %$% rbeta(1000, a, b)
mean(rvec)
sd(rvec)
plot(xs, dens, type='l')

# problem 10
get_posterior_mean <- function(xs, posterior) 
  sapply(xs, function(n) n * posterior(n)) %>% sum()

prior <- function(N) {
  (1/100) * (99/100)^(N - 1)
}
xs <- 1:500
dens <- sapply(xs, prior)
plot(xs, dens, type='l')

number_seen <- 5
p_data_given_n <- function(n) if (n < number_seen) 0 else 1 / n
prior_on_n <- function(n) (1/100) * (99/100)^(n-1)
posterior <- function(n) p_data_given_n(n) * prior_on_n(n)
max_n <- 1000
xs <- 0:max_n
dens <- sapply(xs, posterior)
c <- 1 / sum(dens)
x_range <- number_seen:max_n
posterior_mean <- sapply(x_range, function(n) c * n * posterior(n)) %>% sum()
posterior_mean
posterior_variance <- x_range %>%
  sapply(function(n) ((n - posterior_mean)^2) * (c/n) * (1/100) * (99/100)^n) %>%
  sum()
posterior_sd <- sqrt(posterior_variance)
posterior_sd
plot(xs, dens, type='l')
abline(v=posterior_mean, col='blue')

# 10c
get_posterior_mean <- function(xs, posterior) {
  sapply(xs, function(n) n * posterior(n)) %>% sum()
}

number_seen <- 203
prior_max <- 1000
xs <- 1:prior_max
# uniform prior
prior <- Vectorize(function(n) if (n <= prior_max) 1/prior_max else 0)
# geometric prior
prior <- Vectorize(function(n) (1/100) * (99/100)^(n - 1))
## These two are broken..
# gaussian prior
# prior <- Vectorize(
#   function(n) {
#     dens <- dnorm(n, 500, 300)
#     dens[dens < 0] <- 0
#     function(n) dens[n]
#   }())
# # exponential prior
# prior <- Vectorize(function(n) {
#   dens <- dexp(xs, 0.2)
#   function(n) dens[n]
# })
# plot(xs, prior(xs), type='l')
p_data_given_n <- Vectorize(function(n) 
  if (n >= number_seen) 1/n 
  else 0
)
unnormalized_posterior <- Vectorize(function(n) 
  prior(n) * p_data_given_n(n)
)
c <- sum(unnormalized_posterior(xs))
posterior_fn <- Vectorize(function(n) unnormalized_posterior(n) / c)
posterior_vec <- posterior_fn(xs)
stopifnot(sum(posterior_vec) == 1)
posterior_mean <- get_posterior_mean(xs, posterior_fn)
plot_settings <- list(fn = plot, color = 'black')
plot_settings <- list(fn = lines, color = 'red')
plot_settings %$% fn(xs, posterior_vec, col=color, type='l')
plot_settings %$% abline(v=posterior_mean, col=color)
posterior_mean
plot_settings %$% lines(xs, prior(xs), type='l', col=color)
abline(v=get_posterior_mean(xs, prior))

# OK so woooo the answer to the taxicab problem definitely depends
# a lot on the prior! Which makes sense cuz one datum but still!
# We could make a fun site where you can play with this problem <hug emoji>
# with different priors, and drawing taxis seen, and, and...
# or a blog post! We can use the evolving-distribution gif style
# from the CLT post. One main insight is that you're moving probability
# from < 203 to the LHS of 203+ with the update.
# Does seeing 503, then 203, differ from seeing 503 alone? Probably no?

#####
# 11
y <- c(43, 44, 45, 46.5, 47.5)
prior <- Vectorize(
  function(theta) if (theta >= 0 && theta <= 100) 1/theta else 0
)
dens <- function (y, thetas) {
  sapply(thetas, function(th) prod(dcauchy(y, th, 1)))
}
dtheta <- .01
theta <- seq(0, 100, dtheta)
unnormalized_posterior <- dens(y, theta)
normalized_posterior <- unnormalized_posterior / sum(unnormalized_posterior)
plot(theta, normalized_posterior, type='l', 
     ylim=c(0, 1.1 * max(normalized_posterior)))
sum(theta * normalized_posterior)

?dcauchy

# 11b
posterior_draws <- sample(theta, 1000, replace=TRUE, prob=normalized_posterior)
hist(posterior_draws, breaks=50)

# 11c
# So to get the posterior predictive distribution for a parameter theta... 
# we take many draws of theta from the posterior... and generate a random
# number from the "likelihood"... is what the website said... so does
# that mean the cauchy is the likelihood here? Hmh.
# WRITE THIS DOWN THOOO :) here's getting closer to knowing how to
# draw from predictive distributions.
posterior_predictive <- rcauchy(length(posterior_draws), posterior_draws, 1)
hist(posterior_predictive, breaks=100)

# 13a

compute_gamma_parameters <- function(mu, var) {
  theta <- ((mu + var) / mu) - 1
  k <- mu / theta
  
  beta <- 1 / theta
  alpha <- k
  list(alpha = alpha, beta = beta)
}

data_str <- "1976 24 734 0.19
1977 25 516 0.12
1978 31 754 0.15
1979 31 877 0.16
1980 22 814 0.14
1981 21 362 0.06
1982 26 764 0.13
1983 20 809 0.13
1984 16 223 0.03
1985 22 1066 0.15"
dat <- read.table(text=data_str, sep = ' ', 
                  col.names = c('year', 'accidents', 'deaths', 'death_rate'))
# P(H | data, X) =~ P(data|X) * P(X)
variable <- list(name='accidents', xs=0:200, prior_mean=30, prior_stddev=27.5)
variable <- list(name='deaths', xs=0:4000, prior_mean=1000, prior_stddev=500)
y <- dat[,variable$name]
gamma_params <- variable %$% compute_gamma_parameters(prior_mean, prior_stddev^2)
prior <- gamma_params %$% dgamma(variable$xs, shape=alpha, scale=1/beta)

plot(variable$xs, prior, type = 'l')

# By conjugacy - a gamma prior with a poisson likelihood 
# gives a gamma posterior with parameters given by these sums.
posterior_alpha <- gamma_params %$% {alpha + sum(y)}
posterior_beta <- gamma_params %$% beta + length(y)

posterior <- dgamma(variable$xs, shape=posterior_alpha, scale=1/posterior_beta)
posterior_draws <- sample(variable$xs, 100000, replace=TRUE, prob=posterior)

# here... we wanna find... the probability of each possible new data point,
# given the location parameter.
# maybe it's: draw a location from the posterior, then create a new
# distribution of the same form as the posterior but with a new location,
# then draw a data point from it? And do this many times?
#
# OK! 
# The posterior here is the posterior _over the location parameter_.
# So I need the probability of the predicted point y` for a bunch of y`s 
# generated from distributions with locations drawn from that posterior.
# So I did it right!
pospred_ndraws <- 100000
locations <- sample(variable$xs, pospred_ndraws, replace=TRUE, prob=posterior)
# The mean is gamma-distributed, but the data are poisson distributed.
# So when doing simulated data draws, use the poisson - because we're drawing
# from the data distribution!
posterior_predictive <- sapply(
  locations, 
  function(location) {
    dens <- rpois(1, location)
    dens
  })

# for comparison
draws_from_posterior <- rgamma(
  pospred_ndraws, shape=posterior_alpha, scale=1/posterior_beta)

hist(draws_from_posterior, breaks=80)

posterior[posterior == Inf] <- 0
posterior %<>% {. / sum(.)}
alpha = 0.4
ggplot() + 
  geom_histogram(aes(x=draws_from_posterior), bins=100, fill = 'red', alpha=alpha) +
  geom_histogram(aes(x=posterior_predictive), bins=30, fill = 'blue', alpha=0.2) +
  theme_bw() 

# these two are same!
qgamma(c(0.025, 0.975), shape=posterior_alpha, scale=1/posterior_beta)
quantile(locations, probs = c(0.025, 0.975))

# and the wider predictive distribution:
quantile(posterior_predictive, probs = c(0.025, 0.975))
# Wow... so if you do have a posterior and want to predict from it,
# you need to do this extra thing! Don't accidentally talk about
# location parameters when you mean to talk about distributions!
sum(variable$xs * prior)
sum(variable$xs * posterior)
mean(y)


# 13b
# getting miles_flown algebraically:
# death_rate = (deaths / miles_flown) * 1e8
# miles_flown = (deaths / death_rate) * 1e8
dat %<>% dplyr::mutate(miles_flown = (deaths / death_rate) * 1e8)
variable <- list(name = "rate", xs = seq(0, 0.01, 0.0000001))
prior <- list(alpha = 0, beta = 0)
posterior_alpha <- (prior$alpha + sum(y))
posterior_beta <- (prior$beta + sum(dat$miles_flown))


miles_flown_to_predict <- 8e10
n_posterior_draws <- 10000
posterior_draws <- rgamma(n_posterior_draws, 
                          shape=posterior_alpha, scale=1/posterior_beta)
# Unlike the next problem, here we are plugging in a single value for 
# the explanatory variable - miles flown - and seeing what the prediction
# of deaths for that number of miles flown is. And this number is hiiiigh!
# So we get a distribution with a higher mean.
posterior_predictive <- rpois(1000, posterior_draws * miles_flown_to_predict)
hist(posterior_predictive)
alpha = 0.4
ggplot() + 
  geom_histogram(aes(x=posterior_draws * miles_flown_to_predict), 
                 bins=100, fill = 'red', alpha=alpha) +
  geom_histogram(aes(x=posterior_predictive), bins=30, fill = 'blue', alpha=0.2) +
  theme_bw()

quantile(posterior_predictive, probs=c(0.025, 0.975))


# 20
# a
# Average over all the posteriors for y >= 100.
xs <- 0:2000
prior <- exp()
#guh 

# 21

# Estimate the percentage of the (adult) population in each state (excluding
# Alaska, Hawaii, and the District of Columbia) who label themselves as ‘very liberal’.
remove_names <- c('alaska', 'hawaii', 'washington dc')
xs <- seq(0, 1, 0.0001)
survey_df <- foreign::read.dta('data/pew_research_center_june_elect_wknd_data.dta')
variable <- "very_liberal"
survey_df[,variable] <- survey_df['ideo'] == 'very liberal'
states_responses <- survey_df %>%
  {base::split(., .[,"state"])}
states_responses %<>% {.[!(names(.) %in% remove_names)]}
yeses <- sapply(states_responses, 
                function(dat) sum(dat[,variable], na.rm = TRUE))
noes <- sapply(states_responses, 
               function(dat) sum(!dat[,variable], na.rm = TRUE))
respondents <- yeses + noes

data_df <- tibble::tibble(
  state = names(states_responses),
  respondents = respondents,
  prop = yeses / (yeses + noes))
prior_alpha <- 3
prior_beta <- 20
prior <- dbeta(xs, prior_alpha, prior_beta)
plot(xs, prior, type='l')

posteriors <- Map(
  function(nyes, nno) dbeta(xs, prior_alpha + nyes, prior_beta + nno) / length(xs),
  yeses, noes)
intervals <- Map(
  function(nyes, nno) 
    # bayestestR::hdi(posterior, ci=0.95),
    qbeta(c(0.025, 0.975), prior_alpha + nyes, prior_beta + nno),
    yeses, noes)

means <- sapply(posteriors, function(p) sum(xs * p))
lbs <- sapply(intervals, function(tuple) tuple[1])
ubs <- sapply(intervals, function(tuple) tuple[2])
intervals_df <- tibble::tibble(state = names(intervals),
                               lb = lbs, ub = ubs, mid = means)
ordering <- intervals_df %>% dplyr::arrange(mid) %$% state
intervals_df$state %<>% factor(levels = rev(ordering))

ggplot(intervals_df, aes(y = state)) +
  geom_segment(aes(yend = state, x = lb, xend = ub)) +
  geom_point(aes(x = mid)) +
  theme_bw()

results_df <- read.csv('data/2008ElectionResult.csv')

respondents_df <- survey_df %>% dplyr::select(state, )
states_df <- tibble::tibble(state = tolower(state.name), state_code = state.abb)

results_df$state %<>% tolower()
together_df <- intervals_df %>%
  dplyr::inner_join(results_df, by = 'state') %>%
  dplyr::inner_join(states_df, by = 'state') %>%
  dplyr::inner_join(data_df, by = 'state')

ggplot(together_df, aes(x = vote_Obama_pct)) +
  geom_segment(aes(xend = vote_Obama_pct, y = lb, yend = ub), 
               color='blue', alpha = 0.10) +
  geom_text(aes(y = mid, label = state_code), color = 'blue') +
  geom_text(aes(y = prop, label = state_code), color = 'red') +
  theme_bw() +
  labs(y = "% very liberal") +
  theme(panel.grid = element_blank())



