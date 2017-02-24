#' Plot the model of an isoform pair
#'
#' @param deltaVal A vector of difference between two isoforms of the isoform pair
#' @param dmix A mixture model from doGaussianMixtureModel function
#' @param plot.title Title of the plot
#' @param cex.main.size Font size of the plot
#' @param fit.line An option for plotting the fitting line of distribution
#' @param group.color A vector of colors for the cell groups
#' @param lwd The thickness of lines
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
#' tbreak=round(sqrt(ncol(isoformDataSample)))
#' x0=isoformDataSample[1,]
#' x1=isoformDataSample[2,]
#' model.res=doMixtureModelIsoformPair(x0,x1,tbreak=tbreak)
#' dmix=model.res$dmix.list[[model.res$AIC.min.Id]];
#' deltaVal=x0-x1
#' plotHistModels(deltaVal=deltaVal,dmix=dmix)

plotHistModels <-
function(deltaVal,dmix,plot.title=NULL,cex.main.size=NULL,fit.line=FALSE, group.color=c("darkgreen","darkorange","darkmagenta","darkorchid"),lwd=1){
	nq=length(dmix$p)
	param =dmix$opt$par

	delta = c(param[1:nq])
	beta = param[(nq+1):(2 * nq)]
	sigma = param[(nq*2+1):(3 * nq)]

	sigma=sigma*sigma
	sigma = ifelse(sigma <0.1,0.1,sigma)


	if (is.null(cex.main.size)) cex.main.size=1
	if (is.null(plot.title)) plot.title=""
  
	x<- seq(min(deltaVal),max(deltaVal),len=100)
	hist(deltaVal,prob=T,main=plot.title,cex.main = cex.main.size,xlab=NULL,ylab=NULL)
#	lines(x,dmix$p0.raw*dnorm(x,delta[1],sigma[1]),col="black")
	for (u in 1:length(dmix$p))
		lines(x,dmix$p[u]*dnorm(x,delta[u],sigma[u]),col=group.color[u],lwd=lwd)

	if (fit.line){
		den=dmix$p[1]*dnorm(x,delta[1],sigma[1])
		for (u in 1:(length(dmix$p)-1)){
			den=den+dmix$p[u+1]*dnorm(x,delta[u+1],sigma[u+1])
		}
		points(x,den,col="black",pch=20)
	}
}
