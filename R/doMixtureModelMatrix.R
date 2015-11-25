#' Run mixture models for a data matrix of isoform expression
#'
#' @param isoformLevel.data A data matrix of isoform expression
#' @param txdb A txdb of annotation reference
#' @param usePlot A boolean value to allow plotting data and results, usePlot=FALSE (default)
#' @param useMixModel A boolean value to allow implementing the mixture model, useMixModel=TRUE (default)
#' @param zero.thres The threshold to separate components from null component p0. See \code{\link{doMixtureModelIsoformPair}} function
#' @param nq.num The maximum components of the mixture model. See \code{\link{doMixtureModelIsoformPair}} function
#' @param fn.out File name of RData to export results
#' @param tbreak The number of breaks to create bins. See \code{\link{doMixtureModelIsoformPair}} function
#' @param group.label The group label of cells used for pair plotting
#' @return Lists of detected mixture models, AIC scores, number of components and state of models
#' @importFrom AnnotationDbi select
#' @export
#' @examples
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' library("AnnotationDbi")
#' txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
#' #The txdb also can be imported from a gtf file by using package GenomicFeatures
#' data(isoformDataSample)
#' #preprocessing
#' isoformDataSample=ifelse(isoformDataSample <= 3,0,isoformDataSample)
#' isoformDataSample=isoformDataSample[which(rowSums(isoformDataSample)>0),]
#' #tranform read count dataset to log scale
#' isoformDataSample=logC(isoformDataSample)
#' #now data is ready
#' tbreak=round(sqrt(ncol(isoformDataSample)))
#' res=doMixtureModelMatrix(isoformLevel.data=isoformDataSample,txdb=txdb,tbreak=tbreak)
doMixtureModelMatrix <-
function(isoformLevel.data,txdb=NULL,usePlot=FALSE,useMixModel=TRUE,zero.thres=0.05,nq.num=3,fn.out=NULL,tbreak=NULL,group.label=NULL){
	isoformNames=rownames(isoformLevel.data)
	isoformQueryRes=select(txdb, keys = isoformNames, columns="GENEID", keytype="TXNAME")

	res=duplicated(isoformQueryRes$GENEID)
	isoformGeneList=unique(isoformQueryRes$GENEID[which(res==TRUE)])
	na.ids=which(is.na(isoformGeneList))
	if (length(na.ids)>0) isoformGeneList=isoformGeneList[-na.ids] #remove NA

	dmix.list=NULL;
	AIC.list=NULL;
	nq.list=NULL;		
	ok.model.list=NULL;

	for (i in 1:length(isoformGeneList)){
		isoMatId=which(isoformQueryRes$GENEID==isoformGeneList[i])				
		isoMat=isoformLevel.data[isoMatId,]		

		res=doMixtureModelGene(isoMat,geneName=isoformGeneList[i],usePlot=usePlot,useMixModel=useMixModel,zero.thres=zero.thres,nq.num=nq.num,tbreak=tbreak,group.label=group.label)

		ok.model.list=c(ok.model.list,res$ok.model.list)
		dmix.list=c(dmix.list,res$dmix.list)
		nq.list=c(nq.list,res$nq.list)
		AIC.list=c(AIC.list,res$AIC.list)
		if (usePlot)	readline(prompt="Press [enter] to continue")
	}

	if (!is.null(fn.out)) save(dmix.list, AIC.list, nq.list, ok.model.list,file=fn.out)
	return(list(dmix.list=dmix.list, AIC.list=AIC.list, nq.list=nq.list, ok.model.list=ok.model.list))
}
