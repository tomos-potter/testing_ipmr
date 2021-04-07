# install development version of ipmr

#if(!require('remotes', quietly = TRUE)) {
#  install.packages("remotes")
#}

#remotes::install_github("levisc8/ipmr", build_vignettes = TRUE)

library(ipmr)

my_data_list = list(s_int     = -2.2,   # coefficients(my_surv_mod)[1]
                    s_slope   = 0.25,  # coefficients(my_surv_mod)[2]
                    g_int     = 0.2,   # coefficients(my_grow_mod)[1]
                    g_slope   = 0.99,  # coefficients(my_grow_mod)[2]
                    sd_g      = 0.7,   # sd(resid(my_grow_mod))
                    r_r_int   = 0.003, # coefficients(my_pr_flower_mod)[1]
                    r_r_slope = 0.015, # coefficients(my_pr_flower_mod)[2]
                    r_s_int   = 0.45,   # coefficients(my_seed_mod)[1]
                    r_s_slope = 0.075, # coefficients(my_seed_mod)[2]
                    mu_fd     = 2,     # mean(recruit_data$size_next)
                    sd_fd     = 0.3)   # sd(recruit_data$size_next)

my_simple_ipm <- init_ipm(sim_gen = "simple",
                          di_dd   = "di",
                          det_stoch = "det")


my_simple_ipm <- define_kernel(
  
  proto_ipm = my_simple_ipm,
  
  # Name of the kernel
  
  name      = "P_simple",
  
  # The type of transition it describes (e.g. continuous - continuous, discrete - continuous).
  # These must be specified for all kernels!
  
  family    = "CC",
  
  # The formula for the kernel. We dont need to tack on the "z'/z"s here.  
  
  formula   = s * g,
  
  # A named set of expressions for the vital rates it includes. 
  # note the use of user-specified functions here. Additionally, each 
  # state variable has a stateVariable_1 and stateVariable_2, corresponding to
  # z and z' in the equations above. We don't need to define these variables ourselves,
  # just reference them correctly based on the way we've set up our model on paper.
  
  # Perform the inverse logit transformation to get survival probabilities
  # from your model. plogis from the "stats" package does this for us. 
  
  s         = plogis(s_int + s_slope * dbh_1), 
  
  # The growth model requires a function to compute the mean as a function of dbh.
  # The SD is a constant, so we don't need to define that in ... expression, 
  # just the data_list.
  
  g         = dnorm(dbh_2, mu_g, sd_g),
  mu_g      = g_int + g_slope * dbh_1,
  
  
  # Specify the constant parameters in the model in the data_list. 
  
  data_list = my_data_list,
  states    = list(c('dbh')),
  
  # If you want to correct for eviction, set evict_cor = TRUE and specify an
  # evict_fun. ipmr provides truncated_distributions() to help. This function
  # takes 2 arguments - the type of distribution, and the name of the parameter/
  # vital rate that it acts on.
  
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun    = 'norm',
                                      target = 'g')
) 

my_simple_ipm <- define_kernel(
  proto_ipm = my_simple_ipm,
  name      = 'F_simple',
  formula   = r_r * r_s * r_d,
  family    = 'CC',
  
  # Inverse logit transformation for flowering probability
  # (because we used a logistic regression)
  
  r_r       = plogis(r_r_int + r_r_slope * dbh_1),
  
  # Exponential function for seed progression 
  # (because we used a Poisson)
  
  r_s       = exp(r_s_int + r_s_slope * dbh_1),
  
  # The recruit size distribution has no maternal effect for size,
  # so mu_fd and sd_fd are constants. These get passed in the 
  # data_list
  
  r_d       = dnorm(dbh_2, mu_fd, sd_fd),
  data_list = my_data_list,
  states    = list(c('dbh')),
  
  # Again, we'll correct for eviction in new recruits by
  # truncating the normal distribution.
  
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun    = 'norm',
                                      target = 'r_d')
) 

# Next, we have to define the implementation details for the model. 
# We need to tell ipmr how each kernel is integrated, what state
# it starts on (i.e. z from above), and what state
# it ends on (i.e. z' above). In simple_* models, state_start and state_end will 
# always be the same, because we only have a single continuous state variable. 
# General_* models will be more complicated.

my_simple_ipm <- define_impl(
  proto_ipm = my_simple_ipm,
  make_impl_args_list(
    kernel_names = c("P_simple", "F_simple"),
    int_rule     = rep("midpoint", 2),
    state_start  = rep("dbh", 2),
    state_end    = rep("dbh", 2)
  )
) 

my_simple_ipm <- define_domains(
  proto_ipm = my_simple_ipm,
  dbh = c(0, # the first entry is the lower bound of the domain.
          50, # the second entry is the upper bound of the domain.
          100 # third entry is the number of meshpoints for the domain.
  ) 
) 

# Next, we define the initial state of the population. We must do this because
# ipmr computes everything through simulation, and simulations require a 
# population state.

my_simple_ipm <- define_pop_state(
  proto_ipm = my_simple_ipm,
  n_dbh = runif(100)
)

my_simple_ipm <- make_ipm(proto_ipm = my_simple_ipm,
                          iterations = 200)


lambda_ipmr <- lambda(my_simple_ipm) 
w_ipmr      <- right_ev(my_simple_ipm, iterations = 200)
v_ipmr      <- left_ev(my_simple_ipm, iterations = 200)


# Next example (from ipmr-introduction)

my_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>% 
  define_kernel(
    name      = "P",
    formula   =  s * g,
    family    = "CC",
    s         = 1/(1 + exp(-(s_int + s_slope * dbh_1))),
    g         = dnorm(dbh_2, g_mu, g_sd),
    g_mu      = g_int + g_slope * dbh_1,
    
    data_list = list(
      s_int     = -3,
      s_slope   = 0.5,
      g_int     = 0.1,
      g_slope   = 1.033,
      g_sd      = 2.2
    ),
    states        = list(c('dbh')),
    has_hier_effs = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun   = "norm", 
                                            target =  "g")
  )

my_ipm <- my_ipm %>% 
  define_kernel(
    name      = "F",
    formula   = r_r * r_s * r_d, # Eq 7
    family    = "CC",
    
    # Eq 8. Again, we use the inverse logit transformation to compute pr(flowering)
    r_r       = 1/(1 + exp(-(r_r_int + r_r_slope * dbh_1))),
    
    # Eq 9. We exponentiate this because of the log link in our seed production model
    r_s       = exp(r_s_int + r_s_slope * dbh_1),
    
    # Eq 10. In this case, both the mean and standard deviation are constants
    
    r_d       = dnorm(dbh_2, r_d_mu, r_d_sd),
    data_list = list(
      r_r_int   = -5,   # coef(my_flower_mod)[1]
      r_r_slope = 0.1,   # coef(my_flower_mod)[2]
      r_s_int   = -3,   # coef(my_seed_mod)[1]
      r_s_slope = 0.03,  # coef(my_seed_mod)[2]
      r_d_mu    = 1.2,   # mean(my_recr_data$size_2, na.rm = TRUE)
      r_d_sd    = 0.7    # sd(my_recr_data$size_2, na.rm = TRUE)
    ),
    states        = list(c('dbh')),
    has_hier_effs = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "r_d")
  )


my_ipm <- my_ipm %>%
  
  define_impl(
    
    # Create one list of lists with all the info
    
    list(
      
      # The components of the big list should be the sub-kernels. Each sub-kernel
      # requires 3 pieces of information: the numerical integration rule, the
      # state variable it modifies (state_start), and the state variable it
      # generates (state_end).
      
      P = list(int_rule  = "midpoint",
               state_start = "dbh",
               state_end   = "dbh"),
      F = list(int_rule  = "midpoint",
               state_start = "dbh",
               state_end   = "dbh")
    )
  )

my_ipm <- my_ipm %>%
  define_domains(
    dbh = c(        # the name of the state variable. 
      1,            # the lower bound for the domain. L from Eq 1
      30,           # the upper bound for the domain. U from Eq 1
      200           # the number of mesh points to use for integration 
    )    
  )

my_ipm <- my_ipm %>%
  define_pop_state(n_dbh = rep(1/200, 200))

my_ipm <- make_ipm(my_ipm,
                   iterations = 100)

lambda_ipmr <- lambda(my_ipm)
repro_value <- left_ev(my_ipm)
stable_dist <- right_ev(my_ipm)

#-----------

# ICEPLANT EXAMPLE
data(iceplant_ex)

# growth model

grow_mod <- lm(log_size_next ~ log_size, data = iceplant_ex)
grow_sd  <- sd(resid(grow_mod))


# survival model

surv_mod <- glm(survival ~ log_size, data = iceplant_ex, family = binomial())

# Pr(flowering) model

repr_mod <- glm(repro ~ log_size, data = iceplant_ex, family = binomial())

# Number of flowers per plant model

seed_mod <- glm(flower_n ~ log_size, data = iceplant_ex, family = poisson())

# New recruits have no size(t), but do have size(t + 1)

recr_data <- subset(iceplant_ex, is.na(log_size))

recr_mu  <- mean(recr_data$log_size_next)
recr_sd  <- sd(recr_data$log_size_next)

# This data set doesn't include information on germination and establishment.
# Thus, we'll compute the realized recruitment parameter as the number
# of observed recruits divided by the number of flowers produced in the prior
# year. 

recr_n   <- length(recr_data$log_size_next)

flow_n   <- sum(iceplant_ex$flower_n, na.rm = TRUE)

# Now, we bundle everything into a list as we did before. We can
# pass the model objects directly into the list, and do not need to extract
# coefficients. 

params <- list(recr_mu = recr_mu,
               recr_sd = recr_sd,
               grow_sd = grow_sd,
               surv_mod = surv_mod,
               grow_mod = grow_mod,
               repr_mod = repr_mod,
               seed_mod = seed_mod,
               recr_n   = recr_n,
               flow_n   = flow_n)

# The lower and upper bounds for integration. Adding 20% on either end to minimize
# eviction. The minimum value is negative, so we multiply by 1.2, rather than 0.8.
# If it were positive, we would need to change that to 0.8.

L <- min(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2
U <- max(c(iceplant_ex$log_size, iceplant_ex$log_size_next), na.rm = TRUE) * 1.2

# now run as before, but subbing in "predict" instead of model formula
pred_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det") %>%
  define_kernel(
    name          = "P",
    family        = "CC",
    formula       = s * g,
    
    # Instead of the inverse logit transformation, we use predict() here.
    # We have to be sure that the "newdata" argument of predict is correctly specified.
    # This means matching the names used in the model itself (log_size) to the names
    # we give the domains (sa_1). In this case, we use "sa" (short for surface area) as
    # the new data to generate predictions for.
    
    s             = predict(surv_mod, 
                            newdata = data.frame(log_size = sa_1), 
                            type    = 'response'),
    g_mu          = predict(grow_mod, 
                            newdata = data.frame(log_size = sa_1),
                            type    = 'response'),
    
    # We specify the rest of the kernel the same way.
    
    g             = dnorm(sa_2, g_mu, grow_sd),
    states        = list(c("sa")),
    data_list     = params,
    has_hier_effs = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun   = "norm", 
                                            target =  "g")
  ) %>%
  define_kernel(
    name          = "F",
    family        = "CC",
    formula       = r_p * r_s * r_d * r_r,
    
    # As above, we use predict(model_object). We make sure the names of the "newdata"
    # match the names in the vital rate model formulas, and the values match the 
    # names of domains they use.
    
    r_p           = predict(repr_mod,
                            newdata = data.frame(log_size = sa_1),
                            type    = 'response'),
    r_s           = predict(seed_mod, 
                            newdata = data.frame(log_size = sa_1), 
                            type    = 'response'),
    r_r           = recr_n / flow_n,
    r_d           = dnorm(sa_2, recr_mu, recr_sd),
    states        = list(c("sa")),
    data_list     = params,
    has_hier_effs = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions(fun   = "norm",
                                            target =  "r_d")
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule     = rep('midpoint', 2),
      state_start    = rep("sa", 2),
      state_end      = rep("sa", 2)
    )
  ) %>%
  define_domains(
    sa = c(L,
           U,
           100)
  ) %>%
  define_pop_state(n_sa = runif(100)) %>%
  make_ipm(iterate = TRUE,
           iterations = 200)


lambda(pred_ipm)

R_0 <- function(ipm) {
  
  Fm <- ipm$sub_kernels$F
  Pm <- ipm$sub_kernels$P
  
  I  <- diag(nrow(Pm))
  
  N  <- solve(I - Pm)
  
  R  <- Fm %*% N
  
  out <- Re(eigen(R)$values[1])
  
  return(out)
  
}

gen_T <- function(ipm) {
  
  R0  <- R_0(ipm)
  lam <- lambda(ipm) %>% as.vector()
  
  out <- log(R0) / log(lam)
  
  return(out)
  
}

R_0(pred_ipm)
gen_T(pred_ipm)
