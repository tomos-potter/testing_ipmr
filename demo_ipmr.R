# install development version of ipmr

#if(!require('remotes', quietly = TRUE)) {
#  install.packages("remotes")
#}

#remotes::install_github("levisc8/ipmr", build_vignettes = TRUE)

# Load ipmr and get the parameter values. The data_list argument for define_kernel
# should hold every regression parameter and every constant used in the model.

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
                    mu_rd     = 2,     # mean(recruit_data$size_next)
                    sd_rd     = 0.3)   # sd(recruit_data$size_next)


my_simple_ipm <- init_ipm(sim_gen = "simple", # simple = single continuosy distributed state variable
                          di_dd   = "di", # density independent
                          det_stoch = "det") # deterministic IPM

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
  # so mu_rd and sd_rd are constants. These get passed in the 
  # data_list
  
  r_d       = dnorm(dbh_2, mu_rd, sd_rd),
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

my_simple_ipm <- make_ipm(proto_ipm = my_simple_ipm)


lambda_ipmr <- lambda(my_simple_ipm)
w_ipmr      <- right_ev(my_simple_ipm)
v_ipmr      <- left_ev(my_simple_ipm, iterations = 200)