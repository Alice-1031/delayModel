library(deSolve)
library(plotly)
library(htmlwidgets)

# Define parameters
amplitude <- c(1,0.1)
mu <- 100
sigma_range <- c(10, 300)
tau_range <- c(0, 300)
initial_p <- 0.5
time_max <- 1000
time_steps <- 10000

# Define functions for s(t) and s'(t)
s <- function(t, A, mu, sigma) {
  ifelse (t<=0,0,A * exp(-((t - mu)^2) / sigma))
}

custom_sin <- function(t) {
  ifelse(t < 0,0,(exp(-0.01*(t-50)) * sin((2*pi*(t-50))/200)))
}

# Define the dual DDE system
dual_dde <- function(t, y, parms) {
  A <- parms$A
  mu <- parms$mu
  sigma <- parms$sigma
  tau <- parms$tau
  interaction <- parms$interaction
  hist_on <- parms$hist_on
  
  p <- y[1]
  
  
  S_lag <- s(t-tau,A,mu,sigma)
  
  
  # hist can be either a custom sine function or zero based on hist_on
  if (hist_on) {hist <- custom_sin(t)} else {hist <- 0}
  # Compute dp/dt based on the interaction flag
  if (interaction) {
    dp_dt <- (S_lag + hist + S_lag * hist) * p * (1 - p)

  } else {
    dp_dt <- (S_lag + hist) * p * (1 - p)
  }

  list(c(dp_dt))
}

# Run simulation for varying tau
run_simulation_for_tau <- function(tau) {
  parms <- list(A = amplitude[1], mu = mu, sigma = sigma_range[1], tau = tau,interaction = FALSE,hist_on = FALSE )
  y0 <- c(p = initial_p)
  times <- seq(0, time_max, length.out = time_steps)

  yy <- dede(y = y0, times = times, func = dual_dde, parms = parms, control = list(mxhist = 1e5))
  sim_df <- as.data.frame(yy)
  sim_df$tau <- tau
  sim_df
}

# Run simulation for varying sigma with amplitude set to 0.1
run_simulation_for_sigma <- function(sigma) {
  parms <- list(A = amplitude[2], mu = mu, sigma = sigma, tau = tau_range[1],interaction = FALSE,hist_on = FALSE)
  
  y0 <- c(p = initial_p)
  times <- seq(0, time_max, length.out = time_steps)

  yy <- dede(y = y0, times = times, func = dual_dde, parms = parms, control = list(mxhist = 1e5))
  sim_df <- as.data.frame(yy)
  sim_df$sigma <- sigma
  sim_df
}

# Run simulation for varying tau with hist
run_simulation_hist <- function(tau) {
  parms <- list(A = amplitude[1], mu = mu, sigma = 100, tau = tau,interaction = FALSE,hist_on = TRUE )

  y0 <- c(p = initial_p)
  times <- seq(0, time_max, length.out = time_steps)

  yy <- dede(y = y0, times = times, func = dual_dde, parms = parms, control = list(mxhist = 1e5,atol = 1e-10, rtol = 1e-10))
  sim_df <- as.data.frame(yy)
  sim_df$tau <- tau
  sim_df
}

# Run simulation for varying tau with hist and interaction term
run_simulation_interaction <- function(tau) {
  parms <- list(A = amplitude[1], mu = mu, sigma = 100, tau = tau,interaction = TRUE,hist_on = TRUE )
  
  y0 <- c(p = initial_p)
  times <- seq(0, time_max, length.out = time_steps)
  
  yy <- dede(y = y0, times = times, func = dual_dde, parms = parms, control = list(mxhist = 1e5,atol = 1e-10, rtol = 1e-10))
  sim_df <- as.data.frame(yy)
  sim_df$tau <- tau
  sim_df
}

# Combine results for all tau values
tau_values <- seq(tau_range[1], tau_range[2], by = 5)
sim_data_tau <- do.call(rbind, lapply(tau_values, run_simulation_for_tau))
sim_data_hist <- do.call(rbind, lapply(tau_values, run_simulation_hist))
sim_data_interaction <- do.call(rbind, lapply(tau_values, run_simulation_interaction))

# Combine results for all sigma values
sigma_values <- seq(sigma_range[1], sigma_range[2], by = 5)
sim_data_sigma <- do.call(rbind, lapply(sigma_values, run_simulation_for_sigma))

# X-axis with range and title
axx <- list(
  range = c(0, time_max),
  title = list(text = "𝑡")
  
)

# Z-axis with range and title
axz <- list(
  range = c(0, 1),
  title = list(text = "𝑝(𝑡)")
)

# Y-axis for Tau plot with range and title
axy1 <- list(
  range = c(0, tau_range[2]),
  title = list(text = "𝜏")  # Italic tau (Unicode)
  
)

# Y-axis for Sigma plot with range and title
axy2 <- list(
  range = c(0, sigma_range[2]),
  title = list(text = "𝜎")  # Italic sigma (Unicode)
  
)

# Tau plot 
tau_plot <- plot_ly(sim_data_tau, x = ~time, y = ~tau, z = ~p, type = "scatter3d", mode = "lines",split = ~tau,line = list(color = "grey"), showlegend = FALSE,scene = 'scene1') 
# Sigma plot 
sigma_plot <- plot_ly(sim_data_sigma, x = ~time, y = ~sigma, z = ~p, type = "scatter3d", mode = "lines",split=~sigma,line = list(color = "grey"), showlegend = FALSE,scene = 'scene2')
#hist_plot
hist_plot <- plot_ly(sim_data_hist, x = ~time, y = ~tau, z = ~p, type = "scatter3d", mode = "lines",split = ~tau,line = list( color = "grey"),showlegend = FALSE,scene = 'scene3')
#interaction_plot
interaction_plot <- plot_ly(sim_data_interaction, x = ~time, y = ~tau, z = ~p, type = "scatter3d", mode = "lines",split = ~tau,line = list( color = "grey"),showlegend = FALSE,scene = 'scene4')

#combined
fig <- subplot(tau_plot,sigma_plot,hist_plot,interaction_plot,nrows = 2,margin = 0.01)
fig <- fig %>% layout(
  margin = list(
    l = 5,  # Left margin
    r = 5,  # Right margin
    t = 30,  # Top margin
    b = 5   # Bottom margin
  ),
  title  = "Allele frequency with changing tau vs variance",
  scene  = list(domain =  list(x = c(0, 0.52), y = c(0.5, 1)),xaxis=axx,yaxis=axy1,zaxis=axz,aspectmode='cube',camera = list(eye = list(x = 0.2,y = -2.2,z = 1.6))),     
  scene2 = list(domain = list(x = c(0.5, 1), y = c(0.5, 1)),xaxis=axx,yaxis=axy2,zaxis=axz,aspectmode='cube',camera = list(eye =  list(x = 0.2,y = -2.2,z = 1.6))),    
  scene3 = list(domain = list(x=c(0,0.52),y=c(0,0.53)),xaxis = axx,yaxis = axy1,zaxis = axz,aspectmode = 'cube',camera = list(eye = list(x = 0.2,y = -2.2, z = 1.6))),
  scene4 = list(domain = list(x=c(0.5,1),y=c(0,0.53)),xaxis = axx,yaxis = axy1,zaxis = axz,aspectmode = 'cube',camera = list(eye = list(x = 0.2,y = -2.2, z = 1.6)))
)
fig

# Save as an HTML file and open in the browser
saveWidget(fig, "allele_frequency_plot.html", selfcontained = TRUE)
browseURL("allele_frequency_plot.html")