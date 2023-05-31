# a revised version of the function kyotil::p.adj.perm()
# the revision involved adding 'drop=FALSE' on line 24
p.adj.perm2 <- function (p.unadj, p.perms, alpha = 0.05) 
{
  stopifnot(ncol(p.perms) == length(p.unadj))
  if (is.null(names(p.unadj))) {
    if (is.null(colnames(p.perms))) {
      names(p.unadj) <- 1:length(p.unadj)
    }
    else {
      names(p.unadj) <- colnames(p.perms)
    }
  }
  B = dim(p.perms)[1]
  m = length(p.unadj)
  mode(p.unadj) <- "numeric"
  which.are.NA <- which(is.na(p.unadj))
  p.unadj.order <- order(p.unadj)
  p.unadj <- p.unadj[p.unadj.order]
  p.perms <- p.perms[, p.unadj.order, drop = FALSE]
  len = sum(round(p.unadj, 2) <= alpha)
  p.FWER = rep(NA, length(p.unadj))
  for (j in 1:len) {
    p.FWER[j] = sum((apply(p.perms[, j:m, drop=FALSE], 1, min, na.rm = T) <= 
                       p.unadj[j]))/B
  }
  p.FWER[1:len] = cummax(p.FWER[1:len])
  p.FDR = rep(NA, length(p.unadj))
  for (j in 1:len) {
    R0_by_resample = apply(p.perms <= p.unadj[j], 1, sum, 
                           na.rm = T)
    ER0 = sum(R0_by_resample)/B
    R.ob = j
    R = sum(pmax(R0_by_resample, R.ob))/B
    p.FDR[j] = min(ifelse(R > 0, ER0/R, 0), 1)
  }
  o1 = order(p.unadj[1:len], decreasing = TRUE)
  ro = order(o1)
  p.FDR[1:len] = pmin(1, cummin(p.FDR[1:len][o1]))[ro]
  .rv <- data.frame(mark=names(p.unadj), p.unadj, p.FWER, p.FDR)
  return(.rv)
}
