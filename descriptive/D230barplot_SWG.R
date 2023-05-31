library(ggplot2)

datDir <- "t:/vaccine/p704/analysis/sieve/tables/DVE"
figDir <- "t:/vaccine/p704/analysis/sieve/figures/descriptive"

dat <- read.csv(file.path(datDir, "HVTN704_site_scanning_results_DVE.csv"))
dat <- subset(dat, hxb2Pos==230)

dat1 <- data.frame(x=c("D230 Match", "D230 Mismatch"), 
                   PE=c(dat$VE.vs.Match.estimate, dat$VE.vs.Mismatch.estimate), 
                   LB=c(dat$VE.vs.Match.CI.low, dat$VE.vs.Mismatch.CI.low),
                   UB=c(dat$VE.vs.Match.CI.high, dat$VE.vs.Mismatch.CI.high))


pdf(file.path(figDir, "D230barplot_SWG.pdf"), width=2.7, height=2.7)
ggplot(dat1, aes(x=x, y=PE)) +
  geom_col(fill="gray50", colour="black") +
  geom_errorbar(aes(ymin=LB, ymax=UB), width=0.3, size=1) +
  geom_text(aes(label=paste0(round(PE, digits=0), "%")), vjust=c(-0.5, 1.5), colour="white") +
  geom_hline(yintercept=0) +
  xlab("") +
  ylab("Prevention Efficacy (%)") +
  scale_y_continuous(breaks=seq(-25, 100, by=25)) +
  coord_cartesian(ylim=c(-35, 100))
dev.off()

