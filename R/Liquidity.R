Liquidity <- function(Fit){
  
  mSmoothedProbs = Fit$mu_all
  mEta = Fit$mEta
  
  vEta_smoothed = apply(mSmoothedProbs * mEta, 2, sum)
  vB_smooted = Fit$vY - vEta_smoothed
  # 
  # vNewComers_Prop = vEta_smoothed/Fit$vY
  # 
  # vLiquidity = vB_smooted/vEta_smoothed
  
  return(list(
    "Newcomers" = vEta_smoothed,
    "Splitters" = vB_smooted
  ))
  
}
