fn.solnp <- function(par0, FUN, ...) {
  
  solver.ctr <- list(trace = 0, rho = 1, outer.iter = 400, inner.iter = 800, delta = 1e-07, tol = 1e-08)
  
  optimiser = suppressWarnings(solnp(par0, FUN, ...))
  
  out = list(pars = optimiser$pars,
             value = tail(optimiser$values, 1),
             hessian = optimiser$hessian,
             convergence = optimiser$convergence)
  
  return(out)
  
}

fn.optim <- function(par0, FUN, ...) {
  
  # solver.ctr <- list(trace = 0, abstol = 1e-8, reltol = 1e-8)
  
  optimiser = suppressWarnings(optim(par0, FUN,
                                     method = "BFGS",
                                     # control = solver.ctr,  
                                     hessian = TRUE, ...))
  
  out = list(pars = optimiser$par,
             value = optimiser$value,
             hessian = optimiser$hessian,
             convergence = optimiser$convergence)
  
  return(out)
  
}

fn.DEoptim <- function(par0, FUN, itermax = 200, ...) {
  
  iM = length(par0)
  
  NP = iM * 10
  
  initialpop = matrix(c(par0, replicate(iM - 1, par0 * (1.0 + rnorm(iM) * 0.1))),
                      nrow = NP, ncol = length(par0), byrow = TRUE)
  
  solver.ctr <- DEoptim::DEoptim.control(itermax = itermax, initialpop = initialpop, NP = NP, cluster = list(...)[["cluster"]])
  
  optimiser = suppressWarnings(DEoptim::DEoptim(FUN, 
                                       control = solver.ctr, ...))
  
  pars = optimiser$optim$bestmem
  value = optimiser$optim$bestval
  
  
  names(pars) = list(...)[["namesPw"]]
  mHessian = NULL
  
  out = list(pars = pars,
             value = value,
             hessian = mHessian,
             convergence = optimiser$convergence)
  
  return(out)
  
}


fn.DEoptim_random <- function(par0, FUN, itermax = 200, ...) {
  
  iM = length(par0)
  
  NP = iM * 10
  
  initialpop = NULL
  
  solver.ctr <- DEoptim::DEoptim.control(itermax = itermax, initialpop = initialpop, NP = NP, cluster = list(...)[["cluster"]])
  
  optimiser = suppressWarnings(DEoptim::DEoptim(FUN, 
                                       control = solver.ctr, ...))
  
  pars = optimiser$optim$bestmem
  value = optimiser$optim$bestval
  
  
  names(pars) = list(...)[["namesPw"]]
  mHessian = NULL
  
  out = list(pars = pars,
             value = value,
             hessian = mHessian,
             convergence = optimiser$convergence)
  
  return(out)
  
}

fn.optim_DEoptim <- function(par0, FUN, itermax = 200, ...) {
  
  optimiser1 = fn.DEoptim(par0, FUN, itermax = 20, ...)
  
  par0 = optimiser1$par
  
  optimiser2 = fn.optim(par0, FUN, ...)
  
  return(optimiser2)
  
}
