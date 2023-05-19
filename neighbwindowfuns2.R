##############################################################
##############################################################
## Script title: neighbourhood window functions 
## Repository title: Testing connectivity indicators in R
##
## Author: Jenny A. Hodgson
## 
## Copyright: (c) Jenny A. Hodgson 2020
##  
## Email: jenny.hodgson@liverpool.ac.uk; github@condatis.org.uk
## 
## Preferred attribution: "derived from" Hodgson (2020) 
## 'Testing connectivity indicators in R' www.github.com/condatis/indicators2020, 
## DOI: [quote DOI given on github page after release] 
## 
##  License: Open Government License v.3 (see license file shared in the repository)
## 
## Disclaimer:The Information is licensed 'as is' and the Information Provider
## and/or Licensor excludes all representations, warranties, obligations and 
## liabilities in relation to the Information to the maximum extent permitted by law. 
## 
## IP NOTE: Any IP introduced in this script is background IP introduced
## by Jenny Hodgson, within the meaning of IP clauses in 
## CWI-Defra-CEH Aligned Research Agreement,Subcontract PROJECT NUMBER 07111,
## Inteum Number 7123, between UK Centre for Ecology & Hydrology and
## University of Liverpool
## 
## Purpose of script: Functions for moving window landscape evaluation
##                    adapted for England analysis to cope with fractional 
##                    cell coverage
##
## Warning: although much of this code has potential to be re-used 
## in projects other than the one for which this script was intended,
## it was not written with broad re-usability in mind. Only
## minimal comments are included for guidance 
##############################################################
##############################################################



hanski<- function(raster,dd,resn,conresn=25,cutoff=0.975,returnrasters=TRUE){
	# Note mask is not an option - we rely on fact that raster is NA outside mask
	# resn must be the resolution of the supplied raster
	# conresn must be the same or smaller. connectivity kernel will be created at this resolution then aggregated
	# dd is dispersal distance in same units as the geography, i.e. usually metres
	# a dispersal kernel theoretically stretches out to infinite distance. It will be truncated so it includes cutoff of its theoretical integrated volume

if(resn<conresn){stop("conresn must be the same or smaller than resn")}

if(identical(resn,conresn)){
	mat<- expweightmat(dd=dd,resn=conresn,cutoff=cutoff)

	conr<- focal(x=raster, w=mat,  na.rm=TRUE, pad = TRUE, padValue = NA)

	#only count connectivity of habitat cells
	conh<- conr*raster

	#average connectivity of habitat, where 1 would be found if 100% coverage
	coni<- cellStats(conh,mean)

}else{ #where conresn is smaller than resn


	aggfact<- ((resn/conresn) %/% 2)*2 +1
	conresn<- resn/aggfact #so conresn may be slightly changed if it does not divide resn by an odd number

	mat<- expweightmat(dd=dd,resn=conresn,cutoff=cutoff)

	matdim<- dim(mat)[1]
	subdists<- (1:matdim)-median(1:matdim)

	#extend the matrix with zeros so that it contains a whole, odd number of larger squares, with the peak of the kernel in the centre
	mat2<-raster(mat, xmn=min(subdists)-0.5, xmx=max(subdists)+0.5, ymn=min(subdists)-0.5, ymx=max(subdists)+0.5)
	size2<- (ceiling(matdim/aggfact) %/% 2)*2 +1
	size<- size2*aggfact
	mat2<- extend(mat2,(size-matdim)/2,value=0)
	#aggregate the connectivity matrix
	coarsemat<- aggregate(mat2,fact=aggfact,fun=sum,na.rm=TRUE)

	conr<- focal(x=raster, w=as.matrix(coarsemat),  na.rm=TRUE, pad = TRUE, padValue = NA)

	#only count connectivity of habitat
	conh<- conr*raster

	#average connectivity of habitat, where 1 would be found if 100% coverage
	coni<- cellStats(conh,mean)

}#end else
if(returnrasters){
return(list(conr=conr,conh=conh,coni=coni))
}else{
return(coni)
}
}

expweightmat<- function(dd,resn,cutoff=0.975){
#all distances should be in metres, if to be used with a raster in m
	# a dispersal kernel theoretically stretches out to infinite distance. It will be truncated so it includes cutoff of its theoretical integrated volume

if(cutoff>=1){stop("cutoff must be less than 1")}

	alpha<- 2/dd

	#find correct cutoff distance, as alpha*r

	rem<- 1-cutoff

	try<-seq(-log(rem),pmax(3,-log(rem)*2), -log(0.95))

	propk<- exp(-try)*(try+1)

	rdd<- try[order(abs(propk-rem))[1]]

	#radius in m
	rd<- rdd/alpha

	#radius in no. of cells
	rc<- round(rd/resn)

	mat<- sqrt(outer((0:rc)^2, (0:rc)^2, FUN = "+"))

	fact<- resn*resn* alpha^2/2/pi
	mat<- fact*exp(-mat*resn*alpha)

	mat[mat< (fact*exp(-rdd)) ]<- 0
	#don't count self
	mat[1,1]<-0

	#construct 4 quadrants
	mat2<- cbind( mat[,(rc+1):2],mat)
	mat4<- rbind( mat2[(rc+1):2,],mat2)

	#rescale to sum to 1
	mat4<- mat4/sum(mat4)

	mat4}


ccvc<- function(raster,resn,rescalc=200,pthresh=0.25,returnrasters=TRUE){
	# Note mask is not an option - we rely on fact that raster is NA outside mask

	# Create matrix show in Taylor,Knight and Harfoot 2014, Figure5
	matq<- matrix(c(0,6,3,6,4,2,3,2,1),nrow=3)
	mat2<- cbind( matq[,3:2],matq)
	mat<- rbind( mat2[3:2,],mat2)
	#rescale to sum to 1
	mat<- mat/sum(mat)

	#aggregation of raster
	aggfact<-round(rescalc/resn)
	if(aggfact>1){#just avoids warning message if aggregation not needed
	coarseraster<- aggregate(raster,fact=aggfact,fun=sum,na.rm=TRUE)
		coarseraster<- (coarseraster/aggfact/aggfact) > pthresh
	}else{
	  coarseraster<- raster> pthresh}#this is assuming that raster of propn coverage will be supplied if it is already at correct resolution
	#convert back to pres/absence using a threshold 


	#do moving window
	conr<- focal(x=coarseraster, w=mat,  na.rm=TRUE, pad = TRUE, padValue = NA)

	#only count connectivity of habitat cells
	conh<- conr*coarseraster
	
	#average connectivity of habitat, where 1 would be found if 100% coverage
	coni<- cellStats(conh,mean)

  
if(returnrasters){
return(list(conr=conr,conh=conh,coni=coni)) 

  
}else{
return(c(coni))
}
}	
