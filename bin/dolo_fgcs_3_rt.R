iterations <- 15000
dt <- 200
checkpoint_interval <- 100
control_interval <- 100
mape_threshold <- rep(0.0035, 13)
zero_abs <- 0.0
#mape_threshold[5] <- 1 #Charge
#ctrl_cell_ids <- seq(0, (400*400)/2 - 1, by = 401)
#out_save <- seq(500, iterations, by = 500)
#out_save = c(seq(1, 10), seq(10, 100, by= 10), seq(200, iterations, by=100))


list(
    timesteps = rep(dt, iterations),
    store_result = FALSE,
    #out_save = out_save,
    checkpoint_interval = checkpoint_interval,
    control_interval = control_interval,
    mape_threshold = mape_threshold,
    zero_abs = zero_abs
)