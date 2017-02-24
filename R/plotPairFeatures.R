#' Plot two isoforms by pair lines
#'
#' @param x0 A vector of data points of the first isoforms
#' @param x1 A vector of data points of the second isoforms
#' @param group.label A vector of group label of cells
#' @param plot.title Title of the plot
#' @param group.color A vector of colors for the cell groups
#' @return None
#' @export
#' @examples
#' data(isoformDataSample)
#' #preprocessing
#' isoformDataSample=ifelse(isoformDataSample <= 3,0,isoformDataSample)
#' isoformDataSample=isoformDataSample[which(rowSums(isoformDataSample)>0),]
#' #tranform read count dataset to log scale
#' isoformDataSample=ifelse(isoformDataSample==0,0,log2(isoformDataSample))
#' #now data is ready
#' x0=isoformDataSample[1,]
#' x1=isoformDataSample[2,]
#' plotPairFeatures(x0,x0)
plotPairFeatures <-
function(x0,x1,group.label=NULL,plot.title=NULL,group.color=c("red","blue","green","pink")){
	if (is.null(plot.title)) plot.title=""
	y0=rep(1.0,length(x0))
	y1=rep(0.0,length(x1))
	xDelta=x0-x1
	if (is.null(group.label)){
	    segCol=ifelse(xDelta>0, group.color[2],group.color[1])
	}else{
	  group.label=as.factor(group.label)
	  segCol=group.color[as.integer(group.label)]
	}
	plot(c(1,1),xlim=c(min(c(x0,x1)),max(c(x0,x1))), ylim=c(0,1), type = "n",yaxt="n",ylab="",xlab="log2 of expression",main=plot.title)
	segments(x0, y0, x1 , y1,col=segCol)
}
