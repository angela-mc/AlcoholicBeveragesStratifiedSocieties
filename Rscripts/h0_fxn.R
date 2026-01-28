

################
#### median ####
################

mediqr<-function(x){
  m <- median(x)
  # use median and association iqr (interquartile range/midspread) because median = value separating the higher half from the lower half of a data sample (prob of values below median = probability values above) ;  given we are summarizing as 95% of values >0 ie summarizing probabilistically, median and associated iqr seem better than mean & sd
  miqr<-IQR(x)
  ymin <- m-miqr/1.35
  ymax <- m+miqr/1.35
  
  toret<-list(m, ymin, ymax)
  names(toret) <- c("med", "ymin", "ymax")
  return(toret)
}

mediqr2<-function(x){
  m <- mean(x)
  # use median and association iqr (interquartile range/midspread) because median = value separating the higher half from the lower half of a data sample (prob of values below median = probability values above) ;  
    # given we are summarizing as 95% of values >0 ie summarizing probabilistically, median and associated iqr seem better than mean & sd
  #miqr<-IQR(x)
  #ymin <- m-miqr/1.35
  #ymax <- m+miqr/1.35
  
  lower95 <- quantile(x, probs = c(0.025, 0.975))[1]
  upper95 <- quantile(x, probs = c(0.025, 0.975))[2]
  lower66 <- quantile(x, probs = c(0.17, 0.83))[1]
  upper66 <- quantile(x, probs = c(0.17, 0.83))[2]
  
  toret<-list(m, lower66, upper66, lower95, upper95)
  names(toret) <- c("mean", "lower66", "upper66", "lower95", "upper95")
  return(toret)
}

####################
#### cond effects ##
####################

condeff_plot <- function(ce_data,title="", addl=F){
  
  allrange2<-max(c(ce_data$estimate__, ce_data$lower__, ce_data$upper__))
  tylim<-c(0,allrange2+0.05 )
  tylim<-c(0,0.85)
  a<-1:5
  cols<-c("darkolivegreen4","cornflowerblue","gold","darkorchid","coral1")
  
  #par(mfrow=c(1,1), mar=c(2.5,6,3.5,1))
  plot(a~1,xaxt = 'n', yaxt = 'n', bty='n', xlab="", pch='', ylim=tylim,ylab="Probability", cex.lab=2, main=title, cex.main=2)
  axis(side=1, tck=0.002, lwd.ticks = 2, line = 0, lwd=1.5, cex.axis=1.5,xlim=tylim, at=c(1,2,3,4,5))
  axis(side=2,tck=0.002, lwd.ticks = 2, line = 0, lwd=1.5, cex.axis=1.5, at=c(0.0,0.2,0.4,0.6,0.8))
  contor=1
  for(k in 1:length(cols)){
    points(x=c(k-0.05, k+0.05),y=c(ce_data$estimate__[contor:(contor+1)]),cex=2, col=cols[k], pch=c(19,17))
    segments(x0=k-0.05,x1=k-0.05, y0=ce_data$lower__[contor], y1=ce_data$upper__[contor], col=cols[k],lwd=2)
    segments(x0=k+0.05,x1=k+0.05, y0=ce_data$lower__[contor+1], y1=ce_data$upper__[contor+1], col=cols[k],lwd=2)
    contor<-contor+2
  }
  if(addl==T){
    legend(x=3.5,y=0.8, legend=c("none", "local com", "local + 1","local + 2", "local + 3"), fill=cols,border=cols,title = "Political complexity", bty="n")
    legend(x=3.5,y=0.6, legend=c("alcohol absent", "alcohol present"), pch=c(19,17), bty="n")
  }
  
}


condeff_plot_ordord <- function(ce_data,title=""){
  
  allrange2<-max(c(ce_data$estimate__, ce_data$lower__, ce_data$upper__))
  tylim<-c(0,allrange2+0.05 )
  tylim<-c(0,0.85)
  a<-1:5
  cols<-c("darkolivegreen4","cornflowerblue","gold","darkorchid","coral1")
  
  par(mfrow=c(1,1), mar=c(2.5,6,3.5,1))
  plot(a~1,xaxt = 'n', yaxt = 'n', bty='n', xlab="", pch='', ylim=tylim,ylab="Probability", cex.lab=1.5, main=title, xlim=c(0,6))
  axis(side=1, tck=0.002, lwd.ticks = 2, line = 0, lwd=1.5, cex.axis=1.5,xlim=tylim, at=c(0,1,2,3,4,5))
  axis(side=2,tck=0.002, lwd.ticks = 2, line = 0, lwd=1.5, cex.axis=1.5, at=c(0.0,0.2,0.4,0.6,0.8))
  contor=1
  for(k in 1:length(cols)){
    points(x=c(k-0.20, k-0.10,k, k+0.10, k+0.20),y=c(ce_data$estimate__[contor:(contor+4)]),cex=2, col=cols[k],bg=cols[k], pch=c(19,15,25,24,23))
    segments(x0=k-0.20,x1=k-0.20, y0=ce_data$lower__[contor], y1=ce_data$upper__[contor], col=cols[k],lwd=2)
    segments(x0=k-0.10,x1=k-0.10, y0=ce_data$lower__[contor+1], y1=ce_data$upper__[contor+1], col=cols[k],lwd=2)
    segments(x0=k,x1=k, y0=ce_data$lower__[contor+2], y1=ce_data$upper__[contor+2], col=cols[k],lwd=2)
    segments(x0=k+0.10,x1=k+0.10, y0=ce_data$lower__[contor+3], y1=ce_data$upper__[contor+3], col=cols[k],lwd=2)
    segments(x0=k+0.20,x1=k+0.20, y0=ce_data$lower__[contor+4], y1=ce_data$upper__[contor+4], col=cols[k],lwd=2)
    contor<-contor+5
  }
  legend(x=3.5,y=0.8, legend=c("none", "local com", "local + 1","local + 2", "local + 3"), fill=cols,border=cols,title = "Political complexity", bty="n")
  legend(x=3.5,y=0.6, legend=c("agr 1", "agr 2", "agr 3", "agr 4", "agr 5"), pch=c(19,15,25,24,23), bty="n", border=cols,pt.bg="black")
}



