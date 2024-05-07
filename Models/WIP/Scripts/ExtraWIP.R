
## pop plotting

if (printpop == TRUE){
  # by age
  p2_msirs <- ggplot(data=run.sim.MSIRS) + geom_line(aes(x=time, y=count, color=age)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    #scale_color_manual(values=colz) +
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRS$time), 1)) + ggtitle("MSIRS") +
    xlab("Year")
  p2_msirs + geom_line(aes(x=time, y=tot_pop))
  p2_msirs
}
if (printpop == TRUE){
  # by class
  agg_adjuv <- aggregate(run.sim.MSIRS$count, by=list(run.sim.MSIRS$time, run.sim.MSIRS$class), FUN=sum)
  colnames(agg_adjuv) = c("time", "class", "count")
  
  p3_msirs <- ggplot(data=agg_adjuv) + geom_line(aes(x=time, y=count, color=class)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRS$time), 1)) + ggtitle("MSIRS") + 
    xlab("Year")
  
  p3_msirs
}


# MSIRN
if (printpop == TRUE){
  # by age
  p2_msirn <- ggplot(data=run.sim.MSIRN) + geom_line(aes(x=time, y=count, color=age)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    #scale_color_manual(values=colz) +
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRN$time), 1)) + ggtitle("MSIRN") +
    xlab("Year")
  p2_msirn + geom_line(aes(x=time, y=tot_pop))
  p2_msirn
  
  # by class
  agg_adjuv <- aggregate(run.sim.MSIRN$count, by=list(run.sim.MSIRN$time, run.sim.MSIRN$class), FUN=sum)
  colnames(agg_adjuv) = c("time", "class", "count")
  
  p3_msirn <- ggplot(data=agg_adjuv) + geom_line(aes(x=time, y=count, color=class)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRN$time), 1)) + ggtitle("MSIRN") + 
    xlab("Year")
  
  p3_msirn
}

# MSIRNI
if (printpop == TRUE){
  # by age
  p2_msirni <- ggplot(data=run.sim.MSIRNI) + geom_line(aes(x=time, y=count, color=age)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    #scale_color_manual(values=colz) +
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRNI$time), 1)) + ggtitle("MSIRNI") +
    xlab("Year")
  p2_msirni + geom_line(aes(x=time, y=tot_pop))
  p2_msirni
  
  # by class
  agg_adjuv <- aggregate(run.sim.MSIRNI$count, by=list(run.sim.MSIRNI$time, run.sim.MSIRNI$class), FUN=sum)
  colnames(agg_adjuv) = c("time", "class", "count")
  
  p3_msirni <- ggplot(data=agg_adjuv) + geom_line(aes(x=time, y=count, color=class)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRNI$time), 1)) + ggtitle("MSIRNI") +
    xlab("Year")
  
  p3_msirni
}

