##############################################################
##############################################################
## Script title: Patch-based functions
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
##  Purpose of script:  function to calculate 'Humphrey Crick's rule
##                       of thumb' and couple of other basic patchwise stats, 
##  			              	starting with a raster of habitat		
## 
## Warning: although much of this code has potential to be re-used 
## in projects other than the one for which this script was intended,
## it was not written with broad re-usability in mind. Only
## minimal comments are included for guidance
##############################################################
##############################################################


library(raster)
library(igraph)
library(rgdal)
library(rgeos)

##############function to characterise patches#############

patchisolation<-function(habitat,sradius,rook=FALSE){

#sradius is search radius for finding neighbour, same units as raster CRS
#habitat raster should be either NA or 0 where there is no habitat
#rook option would make fewer patches, and consequently each patch more likely to have a neighbour (because diagonally-adjacent cells not considered part of same patch)
#rook option may be a better approximation to "same patch" if the habitat starts off at a coarser resolution (but no patch id-ing is going to be perfect at a coarse resolution)
 
resn<-mean(res(habitat))

if(rook){
clumps1<-clump(habitat,directions=4)
}else{
clumps1<-clump(habitat,directions=8)
}

patchsize<-as.data.frame(freq(clumps1))

#everything not a patch is NA  -remove the count of these
patchsize<- patchsize[!is.na(patchsize[,1]),]
names(patchsize)[1]<-"id"
patchsize$areakm<-patchsize$count*resn*resn/1000/1000

#look for presence of neighbours as variance of patch id within sradius
fw<-focalWeight(habitat, d=sradius, type='circle')

focal1<-focal(x=clumps1, w=fw, fun=var, na.rm=T, pad=T,padValue=NA)

#cut out parts of patches where they have a neighbour
patchwneigh<- clumps1*(focal1>0)
#list id numbers of patches if they have a neighbour
idswneigh<-unique(patchwneigh, na.last=NA)

patchsize$hasneigh<-patchsize[,1] %in% idswneigh

return(patchsize)#return this patchwise info so it can be used in different ways
}

################convenience summary function##############

thumb1<- function(patchsize){
list(
withneigh=mean(patchsize$hasneigh),#crick rule of thumb - what proportion of patches have a neighbour?
npatch=dim(patchsize)[1],#number of patches
meanparea=mean(patchsize$areakm)#mean patch area
)
}
