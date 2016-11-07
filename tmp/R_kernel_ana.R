t<-read.table("/Users/peterpfleiderer/Documents/0p5_observed/tmp/tmp_p_to_R.dat") 
pdf=density(t$V1,weights=t$V2/sum(t$V2),from=5.82422118187,to=57.3685089111, kernel="gaussian",bw=1.28860719323,na.rm=TRUE) 
out <- data.frame(x=pdf$x,y=pdf$y/sum(pdf$y))
write.table(out,"/Users/peterpfleiderer/Documents/0p5_observed/tmp/tmp_processed_R_to_p.dat" ,row.names = FALSE,col.names = FALSE ) 
