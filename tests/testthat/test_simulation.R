# Simulation study
library(lavaan)
library(plyr)

set.seed(38327)



# Generate data
get_data <- function(n, m, bNMAR) {
  nclus <- n/m
  stopifnot(floor(nclus) == nclus)

  clus <- rep(1:nclus, each=m)
  
  bx0 <- rt(n = nclus, df = 4)
  by0 <- rt(n = nclus, df = 4)
  by1 <- rt(n = nclus, df = 4) + 1
  
  x <- bx0[clus] + rt(n = n, df = 4) 
  y <- by0[clus] + by1[clus] * x
  
  dat <- cbind(x=x, y=y, clus=clus)
  
  # Missingness is nonignorable wrt first order stats
  missing_prob <- plogis(-3 + bNMAR * (x + y))
  mis <- rbinom(n = 2*n, size = 1, prob = missing_prob)
  dat_mis <- dat
  dat_mis[which(mis == 1)] <- NA
  dat_mis <- as.data.frame(dat_mis)
  
  list(dat_mis=dat_mis, dat_orig=as.data.frame(dat))
}

# Estiamte model
get_estimates <- function(dat) {
  fit <- sem("x~~y", data = dat, meanstructure = TRUE, missing = "fiml")
  stopifnot(lavInspect(fit, "converged"))
  
  est <- coef(fit)
  var_iid <- vcov(fit)
  
  ids <- dat$clus
  empty_idx <- fit@Data@Mp[[1]]$empty.idx
  if(length(empty_idx) > 0) ids <- ids[ -empty_idx ]
  
  var_clus <- vcov_complex(fit, ids = ids)

  res <- c(est=est, var_iid=lav_matrix_vech(var_iid), 
           var_clus=lav_matrix_vech(var_clus))
  res
}


nsim <- 1000

conditions <- expand.grid(bNMAR = c(0, 0.3), n = c(100, 500, 1000, 2000), 
                          m = c(2, 5, 20))
  
#Spop <- cov(subset(get_data(n = 1e7, m = 5, bNMAR=0.3)$dat_orig, select = -clus))
est_pop <- coef(sem("x~~y", sample.cov = Spop, sample.nobs = 1e6))


res <- dlply(conditions, .(bNMAR, n, m), .fun = function(cond) {
  print(cond)
  
  res <- ldply(1:nsim, .fun = function(isim) {
    dat_sim <- get_data(cond[,'n'], cond[,'m'], cond[,'bNMAR'])$dat_mis
    get_estimates(dat_sim) 
  }, .progress = "text")
  
  avg_ests <- apply(res[, grep("^est.", colnames(res))], 2, mean)
  sd_ests <- apply(res[, grep("^est.", colnames(res))], 2, sd)
  #med_ests <- apply(res[, grep("^est.", colnames(res))], 2, median)
  #mad_ests <- apply(res[, grep("^est.", colnames(res))], 2, mad)
  
  #med_vars_iid <- apply(res[, grep("^var_iid", colnames(res))], 2, median)
  #med_vars_clu <- apply(res[, grep("^var_clu", colnames(res))], 2, median)
  avg_vars_iid <- apply(res[, grep("^var_iid", colnames(res))], 2, mean)
  avg_vars_clu <- apply(res[, grep("^var_clus[0-9]", colnames(res))], 2, mean)
  
  #med_se_iid <- sqrt(diag(lavaan::lav_matrix_vech_reverse(med_vars_iid)))
  #med_se_clu <- sqrt(diag(lavaan::lav_matrix_vech_reverse(med_vars_clu)))
  
  avg_se_iid <- sqrt(diag(lavaan::lav_matrix_vech_reverse(avg_vars_iid)))
  avg_se_clu <- sqrt(diag(lavaan::lav_matrix_vech_reverse(avg_vars_clu)))
  
  
  data.frame(
    pop=c(est_pop, 0,0), 
    avg_ests,
    sd_ests,
   # med_ests,
  #  mad_ests,
    avg_se_iid,
    avg_se_clu,
  #  med_se_iid,
  #  med_se_clu,
    avg_se_iid/sd_ests,
    avg_se_clu/sd_ests
  #  med_se_iid/mad_ests,
   # med_se_clu/mad_ests
  )

})


