read_rpkmfiles <- function (infiles) {
    nF<-length(infiles)
    # read data from each file:
    R<-mat.or.vec(0,0)
    C<-mat.or.vec(0,0)
    samples<-mat.or.vec(0,0)
    normreads<-mat.or.vec(0,0)
    allreads<-mat.or.vec(0,0)
    for (i in 1:nF) { 
    	rl<-readLines(con = infiles[i], n = 3) 
    	s <- strsplit(rl[1],"\t")  
    	nS<-length(s[[1]])-1
    	sampl <- s[[1]][2:(nS+1)]
    	s <- strsplit(rl[2],"\t")
    	allr<-as.numeric(s[[1]][2:(nS+1)])
    	s <- strsplit(rl[3],"\t")
    	normr<-as.numeric(s[[1]][2:(nS+1)])
    	samples <- append(samples,sampl)
    	normreads<-append(normreads,normr)
    	allreads<-append(allreads,allr)

    	# read Rpkm
    	data<-read.table(infiles[i], as.is=T, header=F,fill=T)
    	d<-data[,3:(nS+2)] 
    	c<-data[,(nS+3):(nS*2+2)]
    	if (ncol(R)==0) {
       	   R<-d
       	   C<-c
        }else {    
           R<-cbind(R,d)
           C<-cbind(C,c)
        }
    }
    genenames<-data[,1]
    txnames<-data[,2]

    list(rpkm=as.matrix(R),counts=as.matrix(C),normreads=normreads,tx=txnames,genes=genenames,samples=samples,allreads=allreads)
}
