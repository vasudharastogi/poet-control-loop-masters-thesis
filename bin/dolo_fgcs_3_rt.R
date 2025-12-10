iterations <- 10000
dt <- 200
chkpt_interval <- 100
stab_interval <- 100
mape_threshold <- rep(0.01, 13)
mape_threshold[5] <- 1 #Charge
zero_abs <- 1e-13
rb_limit <- 3
rb_interval_limit <- 300
ctrl_cell_ids <- seq(0, (400*400)/2 - 1, by = 401)
out_save <- seq(500, iterations, by = 500)
#out_save = c(seq(1, 10), seq(10, 100, by= 10), seq(200, iterations, by=100))


list(
    timesteps = rep(dt, iterations),
    store_result = TRUE,
    out_save = out_save,
    chkpt_interval = chkpt_interval,
    stab_interval = stab_interval,
    mape_threshold = mape_threshold,
    zero_abs = zero_abs,
    rb_limit = rb_limit,
    rb_interval_limit = rb_interval_limit,
    ctrl_cell_ids = ctrl_cell_ids
)
