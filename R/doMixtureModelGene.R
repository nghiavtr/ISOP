#' Do mixture modelling for all pairs of isoforms of the same single individual gene
#'
#' @param isoMat Data matrix contains isoforms of a same individual gene. Each row indicates a isoform, while each column presents a cell
#' @param geneName A string contains the name of gene that the isoforms are sharing
#' @param usePlot A boolean value to allow plotting data and results. By default usePlot=FALSE 
#' @param useMixModel A boolean value to allow implementing the mixture model. By default useMixModel=TRUE
#' @param zero.thres The threshold to separate components from null component p0. See \code{\link{doMixtureModelIsoformPair}} function
#' @param nq.num The maximum components of the mixture model. See \code{\link{doMixtureModelIsoformPair}} function
#' @param tbreak The number of breaks to create bins. See \code{\link{doMixtureModelIsoformPair}} function
#' @param group.label The group label of cells used for pair plotting
#' @return Lists of results of mixture models for isoform paires
#' @export
#' @examples
#' data(isoformDataSample)
#' #preprocessing
#' isoformDataSample=ifelse(isoformDataSample <= 3,0,isoformDataSample)
#' isoformDataSample=isoformDataSample[which(rowSums(isoformDataSample)>0),]
#' #tranform read count dataset to log scale
#' isoformDataSample=logC(isoformDataSample)
#' #now data is ready
#' tbreak=round(sqrt(ncol(isoformDataSample)))
#' isoMat=isoformDataSample[c(1,2),]
#' res=doMixtureModelGene(isoMat,geneName="118424",tbreak=tbreak)
doMixtureModelGene <-
function(isoMat,geneName="geneName",usePlot=FALSE,useMixModel=TRUE,zero.thres=0.05,nq.num=3,tbreak=NULL,group.label=NULL){
	dmix.list=list();
	AIC.list=list();
	nq.list=list();
	ok.model.list=list();

    if (usePlot)	par(mfrow=c(nrow(isoMat),nrow(isoMat)),mar=c(1,2,2,1)+0.2,oma=c(1,1,1,1));#,cex=1.3)  	
  	for (j in 1:nrow(isoMat)){
  		for (k in 1:nrow(isoMat))
  		if (j!=k){	  			
  			x0=isoMat[j,]
  			x1=isoMat[k,]
  			deltaVal=x0-x1;

  			isOverlapp=FALSE
  			if (sum(abs(deltaVal))==0) isOverlapp=TRUE

  			if (j<k){ 
  				if (usePlot){
  					if (isOverlapp){ plot(c(1,1)); next}
  					plotPairFeatures(x0,x1,group.label=group.label)
  				}
  			}else{ 
  				if (isOverlapp){ 
  					if (usePlot)	plot(c(1,1)); 
  					next
  				}		

  				if (usePlot & !useMixModel) hist(deltaVal,prob=T,main="")

  				if (useMixModel) {
				
  					#run mixture model
  					model.res=doMixtureModelIsoformPair(x0,x1,nq.num=nq.num,tbreak=tbreak,zero.thres=zero.thres)

  					aName=paste(geneName,"^",rownames(isoMat)[j],"^",rownames(isoMat)[k],sep="")
					ok.model.list[[aName]]=model.res$isOk;
					dmix.list[[aName]]=model.res$dmix.list;
					nq.list[[aName]]=model.res$AIC.min.Id
					AIC.list[[aName]]=model.res$AIC.score[model.res$AIC.min.Id]
					if (usePlot){
						dmix=model.res$dmix.list[[model.res$AIC.min.Id]];
						plotHistModels(deltaVal,dmix,plot.title=paste(round(model.res$AIC.score[model.res$AIC.min.Id],2)," (",model.res$AIC.min.Id,")",sep=""),cex.main.size=4/nrow(isoMat),fit.line=FALSE)
					}

				}# end of useMixModel
			}
  		}else{
  			if (usePlot){
	  			plot(c(10,10),type="n",axes=FALSE,frame.plot=TRUE,xlab="",ylab="")  
	  			legendText=rownames(isoMat)[j];
	  			legendCol="black"
	  			legend("center",legendText,text.col=legendCol,cex = 4/nrow(isoMat),bty = "n")
	  		}
  		}
  	}	

  	res=list(dmix.list=dmix.list,AIC.list=AIC.list,nq.list=nq.list,ok.model.list=ok.model.list)
  	return (res)
  }
