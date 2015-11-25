#' @title Get log scale of read counr vector or matrix
#'
#' @description Return log2, log10 or log values of expressed values while assign specific value (log.zero) to zero values
#'
#' @param x A vector or matrix of read counts
#' @param type Type of transformation log2, log10 or log
#' @param log.zero Value that zero values of read counts will become in log scale (=0 by default)
#' @return Data values after log scaling
#' @export
#' @examples
#' data(isoformDataSample)
#' #preprocessing
#' isoformDataSample=ifelse(isoformDataSample <= 3,0,isoformDataSample)
#' isoformDataSample=isoformDataSample[which(rowSums(isoformDataSample)>0),]
#' #tranform read count dataset to log scale
#' isoformDataSample=logC(isoformDataSample)
#' #now data is ready

logC<-function(x,type="log2",log.zero=0){
	#if x is a vector
	if (is.vector(x)){
		xi=x;
		nonZeroIds=which(xi!=0);
		if (length(nonZeroIds)==0) {
			 x[]=log.zero
		} else{
			if (type=="log2") x[nonZeroIds]=log2(x[nonZeroIds])
			if (type=="log10") x[nonZeroIds]=log10(x[nonZeroIds])
			if (type=="log") x[nonZeroIds]=log(x[nonZeroIds])
			x[-nonZeroIds]=log.zero;
		}
		return(x);
	}
	#if x is a matrix
	for (i in 1:ncol(x)){
		xi=x[,i];
		nonZeroIds=which(xi!=0);
		if (length(nonZeroIds)==0) {
			x[,i]=log.zero; 
		} else{
			if (type=="log2") x[nonZeroIds,i]=log2(x[nonZeroIds,i])
			if (type=="log10") x[nonZeroIds,i]=log10(x[nonZeroIds,i])
			if (type=="log") x[nonZeroIds,i]=log(x[nonZeroIds,i])
			x[-nonZeroIds,i]=log.zero;
		}
	}
	return(x)
}