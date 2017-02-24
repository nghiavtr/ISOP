#' Cluster cells to components based on their proximities to the distributions
#'
#' @param iso.pair.names A vector of names of isoform pairs in the format: gene.isoform1.isoform2
#' @param dmix.list A list of mixture models obtained from doMixtureModelMatrix function
#' @param nq.list A list of number of components of detected models
#' @param isoformLevel.data A data matrix of isoform expression
#' @param geneLevel.data A data matrix of gene expression
#' @return Two lists of matrices cellCompMat.conf and cellCompMat.prob indicating confidence scores and probabilities of cells assigned to components
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
#' isoformDataSample=ifelse(isoformDataSample==0,0,log2(isoformDataSample))
#' #now data is ready
#' tbreak=round(sqrt(ncol(isoformDataSample)))
#' model.res=doMixtureModelMatrix(isoformLevel.data=isoformDataSample,txdb=txdb,tbreak=tbreak)
#' iso.pair.names=names(model.res$nq.list)
#' res=assignCellComp(iso.pair.names,model.res$dmix.list,
#'     model.res$nq.list, isoformLevel.data=isoformDataSample)
assignCellComp=function(iso.pair.names,dmix.list,nq.list,isoformLevel.data,geneLevel.data=NULL){	
	cellCompMat.conf=NULL;
	cellCompMat.prob=NULL;
	res=lapply(iso.pair.names, function(x){
		myisopairname=x
		mydmixRes=dmix.list[[myisopairname]]
		mynq=nq.list[[myisopairname]]
		mydmix=mydmixRes[[mynq]]
		s.names=unlist(strsplit(myisopairname,"\\^"))
		if (length(s.names)==3){
			g.name=s.names[1]	
			iso1.name=s.names[2]
			iso2.name=s.names[3]
		} else{
			cat("\n Warning: weird name ",myisopairname," should recheck ... ")
			g.name=paste(s.names[1],"^",s.names[2],sep="")
			iso1.name=s.names[3]
			iso2.name=s.names[4]
		}
		if (g.name != iso2.name){
		  x0= isoformLevel.data[which(rownames(isoformLevel.data)==iso1.name),]
		  x1= isoformLevel.data[which(rownames(isoformLevel.data)==iso2.name),]
		}else{
		  x0= isoformLevel.data[which(rownames(isoformLevel.data)==iso1.name),]
		  x1= geneLevel.data[which(rownames(geneLevel.data)==iso2.name),]
		}
		deltaVal=x0-x1;

		nq=length(mydmix$p)
		param =mydmix$opt$par

		delta = c(param[1:nq])
		beta = param[(nq+1):(2 * nq)]
		sigma = param[(nq*2+1):(3 * nq)]
		sigma=sigma*sigma
		sigma = ifelse(sigma <0.1,0.1,sigma)

		mytbreak=mydmix$tbreak
		mytbreak[1]=min(deltaVal)-1e-6
		mytbreak[length(mytbreak)]=max(deltaVal)+1e-6
		y=cut(deltaVal, mytbreak)
		labs=levels(y)
		breakPoints=cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ), upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))

		assigned.conf=matrix(0,nrow=length(delta),ncol=length(deltaVal),byrow = TRUE)
		colnames(assigned.conf)=names(deltaVal)
		rownames(assigned.conf)=paste(myisopairname,c(1:length(delta)),sep=".")

		for (i in 1:nrow(breakPoints)){
			low.conf=pnorm(breakPoints[i,1],delta,sigma)-0.5
			high.conf=pnorm(breakPoints[i,2],delta,sigma)-0.5
			interval.conf=abs(low.conf)
			#select the closer
			interval.conf=ifelse(interval.conf > abs(high.conf), abs(high.conf), interval.conf)
			#if the component stays inside the interval, all cells in the interval are assigned to the components
			interval.conf[which(sign(low.conf)*sign(high.conf)<0)]=0
			#update to matrix
			myids=which(y==labs[i])
			if (length(myids)>0) assigned.conf[,myids]=interval.conf
		}		
		assigned.conf
	})
	cellCompMat.conf=res;
	names(cellCompMat.conf)=iso.pair.names

	res=lapply(cellCompMat.conf, function(x){
		assigned.conf=x
		assigned.prob=1-assigned.conf
		cSum.tmp=colSums(assigned.prob)
		assigned.prob=apply(assigned.prob,1,function(x) x/cSum.tmp)# already transposed
		assigned.prob
	})
	cellCompMat.prob=res
	
	return(list(cellCompMat.conf=cellCompMat.conf,cellCompMat.prob=cellCompMat.prob))

}