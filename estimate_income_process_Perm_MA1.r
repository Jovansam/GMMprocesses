# The model we simulate and estimate here is an Permanent-transitory model where the permanent component is a random walk and the transitory component is a MA(1)
# Formally: y_it = u_{it} + v_{it} where:
  # u_{it} = u_{it-1} + w_{it}
  # v_{it} = e_{it-1} + \theta e_{it}


# Simulate data
# -------------

rm(list=ls())
set.seed(337631)

library(gmm)
library(data.table)

N <- 10000
T <- 5

sigw <- 0.014
sigu1 <- 0.866
sige <- 0.077
theta <- 0.498

dt <- data.table(n = rep(1:N, each=T), t = rep(1:T, N))
setkey(dt, n, t)


# Random walk component
dt[, w := ifelse(t == 1, 0, rnorm(N*T, mean = 0, sd = sigw))]
dt[, u := ifelse(t == 1, rnorm(N*T, mean = 0, sd = sigu1), 0)]
dt[, u := cumsum(u + w), by=n]

# Moving average component
dt[, e0 := ifelse(t == 1, rnorm(N*T, mean = 0, sd = sige), 0)]
dt[, e := rnorm(N*T, mean = 0, sd = sige)]
dt[, v := e + (theta*shift(e)), by=n]
dt[t == 1, v := e + (theta*e0)]

# Outcome
dt[, y := u + v]


# Moment evaluators
# -----------------

# Function to calculate moments in first difference specification
gfirstdiff <- function(params, x) {
  
  xdt <- as.data.table(x)
  
  tmin <- xdt[, min(t)]
  tmax <- xdt[, max(t)]
  numperiods <- tmax - tmin + 1
  numindivs <- nrow(xdt) / numperiods   # Assumes a balanced panel
  nummoments <- (numperiods - 1) + (numperiods - 2) + (numperiods - 3)
  
  if (numperiods < 4) {
    stop("Too few periods (minimum = 4)")
  }
  
  
  # Theoretical population values
  # -----------------------------
  
  # Order of parameters: sigw, sige, theta
  
  pop_var_dy <- (params[1]^2) + (2.0 * ((params[3]^2) - params[3] + 1) * (params[2]^2))
  pop_cov_dy_Ldy <- -((params[3] - 1)^2) * (params[2]^2)
  pop_cov_dy_L2dy <- -params[3] * (params[2]^2)
  
  
  # Sample estimates
  # ----------------
  
  # First difference outcome
  xdt[, dy := y - shift(y), by=n]
  # De-mean first difference
  xdt[, dy := dy - mean(dy, na.rm = TRUE), by=t]
  
  # Variance
  xdt[, var_dy := dy*dy]

  # Covariance with lagged value
  xdt[, cov_dy_Ldy := dy * shift(dy), by=n]
  
  # Covariance with value lagged two periods
  xdt[, cov_dy_L2dy := dy * shift(dy, n=2), by=n]
  
  # Create matrix to hold moments
  moments <- matrix(NA, nrow = numindivs, ncol = nummoments)

  
  # Calculate moments
  # -----------------
  
  ixmoment <- 1
  for (period in (tmin+1):tmax) {
    
    # Variance
    moments[, ixmoment] <- xdt[t == period, var_dy] - pop_var_dy
    ixmoment <- ixmoment + 1
    
    # Covariance with lagged value
    if (period > tmin + 1) {
      moments[, ixmoment] <- xdt[t == period, cov_dy_Ldy] - pop_cov_dy_Ldy
      ixmoment <- ixmoment + 1
    }

    # Covariance with value lagged two periods
    if (period > tmin + 2) {
      moments[, ixmoment] <- xdt[t == period, cov_dy_L2dy] - pop_cov_dy_L2dy
      ixmoment <- ixmoment + 1
    }
    
  }
  print(params)
  print(colSums(moments))
  return(moments)
  
}



# Function to calculate moments in levels specification
glevels <- function(params, x) {
  
  xdt <- as.data.table(x)
  
  tmin <- xdt[, min(t)]
  tmax <- xdt[, max(t)]
  numperiods <- tmax - tmin + 1
  numindivs <- nrow(xdt) / numperiods   # Assumes a balanced panel
  nummoments <- 0.5 * numperiods * (numperiods + 1)
  
  if (numperiods < 3) {
    stop("Too few periods (minimum = 3)")
  }
  
  # De-mean outcome
  xdt[, y := y - mean(y, na.rm = TRUE), by=t]
  
  # Create matrix to hold moments
  moments <- matrix(NA, nrow = numindivs, ncol = nummoments)
  
  
  # Parameters: sigw, sigu1, sige, theta
  
  # Calculate non-linear combinations of parameters
  pop_var_v <- (1.0 + (params[4]^2)) * (params[3]^2)
	pop_cov_v_Lv = params[4] * (params[3]^2)
  
	ixmoment <- 1
	
  # Loop across diff (between t and s)
	for (diff in 0:(numperiods-1)) {
	
	  # Calculate population covariance of v for current diff
	  pop_cov_v_diff <- 0.0
	  if (diff == 0) pop_cov_v_diff <- pop_var_v
		if (diff == 1) pop_cov_v_diff <- pop_cov_v_Lv
		
		# Calculate sample covariance for diff
		xdt[, cov_y_diff := y * shift(y, n=diff), by=n]
		

		# Loop across period
		for (period in (tmin+diff):tmax) {
		  
		  # Calculate population covariance of u for current t and diff
		  pop_cov_u_t_diff = (params[2]^2) + ((period - diff - tmin) * (params[1]^2))
		  
		  moments[, ixmoment] <- xdt[t == period, cov_y_diff] - pop_cov_v_diff - pop_cov_u_t_diff

		  ixmoment <- ixmoment + 1
		  
		}
				
	}
	
	print(params)
	print(colSums(moments))
	return(moments)
	
}


# Estimate models
# ---------------

reslev <- gmm(glevels, as.matrix(dt), c(sigw = 0.1, sigu1 = 0.5, sige = 0.1, theta = 0.7),
              optfct="nlminb", lower=c(0.001, 0.001, 0.001, 0.1), upper=c(0.4, 0.99, 0.4, 0.99),
              vcov = "iid", prewhite = FALSE, traceIter = TRUE)

resfd <- gmm(gfirstdiff, as.matrix(dt), c(sigw = 0.1, sige = 0.1, theta = 0.7),
             optfct="nlminb", lower=c(0.001, 0.001, 0.1), upper=c(0.4, 0.4, 0.99),
             vcov = "iid", prewhite = FALSE, traceIter = TRUE)

summary(reslev)
confint(reslev,level=.95)

summary(resfd)
confint(resfd,level=.95)


