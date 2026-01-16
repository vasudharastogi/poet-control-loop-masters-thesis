iterations <- 10000
dt <- 200
ctrl_interval <- 200
mape_threshold <- rep(7e-3, 13)
mape_threshold[5] <- 1 #Charge
zero_abs <- 1e-13
rb_limit <- 0
rb_aging_limit <- 2*ctrl_interval

#ctrl_cell_ids <- seq(0, (400*400)/2 - 1, by = 401)
out_save <- c(seq(1, 10), seq(10, 100, by= 10), seq(500, iterations, by = 500))
#out_save = c(seq(1, 10), seq(10, 100, by= 10), seq(200, iterations, by=100))


list(
    timesteps = rep(dt, iterations),
    store_result = TRUE,
    out_save = out_save,
    chkpt_interval = ctrl_interval,
    ctrl_interval = ctrl_interval,
    mape_threshold = mape_threshold,
    zero_abs = zero_abs,
    rb_limit = rb_limit,
    rb_aging_limit = rb_aging_limit
)