#' Do mixture modelling for two isoforms
#'
#' @param x0 A vector of data points from isoform 1
#' @param x1 A vector of data points from isoform 2
#' @param nq.num The maximum components of the mixture model
#' @param tbreak The number of breaks to create bins
#' @param zero.thres The threshold to separate components from null component p0
#' @param min.prop Minimum proportion threshold for a component in the mixture model
#' @return A list consisting of list of mixture models (dmix.list), index of the best model (AIC.min.Id), vector of AIC scores (AIC.score), ...
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
#' res=doMixtureModelIsoformPair(x0,x1,tbreak=tbreak)
doMixtureModelIsoformPair <-
function(x0,x1,nq.num=3,tbreak=NULL,zero.thres=0.05,min.prop=0.025){
	deltaVal=x0-x1;

	deltaProp=c(sum(deltaVal<0),sum(deltaVal==0),sum(deltaVal>0))/length(deltaVal)
	exp.prob=1 - deltaProp[2]
	dmix.list=list();  				
	AIC.score=NULL;
	for (nq in 1:nq.num){
		s0 = exp.prob
		if (nq==1){ 
			s=sd(deltaVal);
			delta=mean(deltaVal);
			res=doGaussianMixtureModel(deltaVal,s=s,delta=delta,useExt=FALSE,tbreak=tbreak)
		}		
		if (nq==2){ 
		# case1: if data are in one side only
			s=c(1,1);
			delta=c(-1,1)			
			quanVal=quantile(deltaVal)
			delta[1]=quanVal[2]
			delta[2]=quanVal[4]
			s[1]=sqrt(abs(quanVal[2]-quanVal[3])/2)
			s[2]=sqrt(abs(quanVal[4]-quanVal[3])/2)
			s=sqrt(s)
			s=ifelse(is.na(s),1,s)
			s=ifelse(s < 0.01,1,s)
			res1=doGaussianMixtureModel(deltaVal,s=s,delta=delta,useExt=FALSE,tbreak=tbreak)

		# case2: if data are in two sides
			s=c(1,1);
			delta=c(-1,1)			
			s[1]=suppressWarnings(sqrt(diff(range(deltaVal[deltaVal < -zero.thres]))/4))
			s[2]=suppressWarnings(sqrt(diff(range(deltaVal[deltaVal > zero.thres]))/4))

			s=ifelse(is.na(s),1,s)
			s=ifelse(s < 0.01,1,s)

			delta[1]=mean(deltaVal[deltaVal< -zero.thres])
			delta[2]=mean(deltaVal[deltaVal>zero.thres])
			delta=ifelse(is.na(delta),0,delta)
			res2=doGaussianMixtureModel(deltaVal,s=s,delta=delta,useExt=FALSE,tbreak=tbreak)
			
			res=res1;
			if (res1$AIC > res2$AIC) res=res2;

		}

		if (nq==3){						
		# case1: if data are in one side only
			s=c(1,1,1);
			delta=c(-1,0,1)								
			quanVal=quantile(deltaVal)
			delta[1]=quanVal[2]
			delta[2]=quanVal[4]
			s[1]=sqrt(abs(quanVal[2]-quanVal[3])/2)
			s[2]=sqrt(abs(quanVal[4]-quanVal[3])/2)
			s=sqrt(s)
			s=ifelse(is.na(s),1,s)
			s=ifelse(s < 0.01,1,s)
			s[3]=exp.prob
			delta[3]=0;
			res1=doGaussianMixtureModel(deltaVal,s=s,delta=delta,useExt=FALSE,tbreak=tbreak)
		# case2: if data are in two sides
			s=c(1,1,1);
			delta=c(-1,0,1)								

			s[1]=suppressWarnings(sqrt(diff(range(deltaVal[deltaVal < -zero.thres]))/4))
			s[2]=suppressWarnings(sqrt(diff(range(deltaVal[deltaVal > zero.thres]))/4))
			s=ifelse(is.na(s),1,s)
			s=ifelse(s < 0.01,1,s)
			delta[1]=mean(deltaVal[deltaVal< -zero.thres])
			delta[2]=mean(deltaVal[deltaVal>zero.thres])
			delta=ifelse(is.na(delta),0,delta)
			s[3]=exp.prob
			delta[3]=0;
			res2=doGaussianMixtureModel(deltaVal,s=s,delta=delta,useExt=FALSE,tbreak=tbreak)
						
			if (deltaProp[1] > 2*min.prop & deltaProp[3]>2*min.prop) res =res2 else{
				if (deltaProp[1] < min.prop | deltaProp[3]<min.prop) res=res1 else {
				res=res1;
				if (res1$AIC > res2$AIC) res=res2;
				}
			}

		}

		if (nq==4){
			s=c(1,1,1,1);
			delta=c(-1,0,0,1)
			quanVal=quantile(deltaVal)
			delta[1]=quanVal[1]
			delta[2]=quanVal[2]
			delta[3]=quanVal[4]
			delta[4]=quanVal[5]
			res=doGaussianMixtureModel(deltaVal,s=s,delta=delta,useExt=FALSE,tbreak=tbreak)
		}

		if (nq>4){
			s=rep(1,nq);
			delta=rep(0,nq);
			tmp=diff(range(deltaVal))/nq
			delta=unlist(lapply(c(1:nq), function(u) tmp*u ))
			res=doGaussianMixtureModel(deltaVal,s=s,delta=delta,useExt=FALSE,tbreak=tbreak)
		}
		dmix.list[[nq]]=res
		AIC.score=c(AIC.score,res$AIC)
	}
	# model selection
	AIC.score.order=order(AIC.score)
	isOk=0;
	AIC.min.Id=AIC.score.order[1]

	for (u in 1:length(AIC.score.order)){
		AIC.min.Id=AIC.score.order[u];
		if (AIC.min.Id==3 & (dmix.list[[AIC.min.Id]]$p0.est!=dmix.list[[AIC.min.Id]]$p0.raw)) next;

		if (sum(dmix.list[[AIC.min.Id]]$p < min.prop)==0){ 
			isOk=1
			break;
		}
	}

	res=list(AIC.score=AIC.score,AIC.min.Id=AIC.min.Id,dmix.list=dmix.list,isOk=isOk)
	return(res)
}
