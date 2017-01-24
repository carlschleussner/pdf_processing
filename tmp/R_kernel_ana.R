t<-read.table("./tmp/tmp_p_to_R.dat") 
pdf=density(t$V1,weights=t$V2/sum(t$V2),from=13.6540403366,to=45.3022842407, kernel="gaussian",bw=0.791206097603,na.rm=TRUE) 
out <- data.frame(x=pdf$x,y=pdf$y/sum(pdf$y))
write.table(out,"./tmp/tmp_processed_R_to_p.dat" ,row.names = FALSE,col.names = FALSE ) 
