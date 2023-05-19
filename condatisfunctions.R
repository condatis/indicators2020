##############################################################
##############################################################
## Script title: condatis functions
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
##  Purpose of script: calculating Condatis speed (optionally flow) 
##                      plus a function to make a study area boundary into  
##  			              appropriate sources and targets			
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


st4directions<-function(loregionb){
  #a function to divide a rasterised study area boundary into 8 segments, which can
  #then be used as sources/targets in Condatis.
  #The input is assumed to be a raster defining a boundary (1s) around a non-disjoint region (0s),
  #surrounded by NAs, for example as produced by the function raster::boundaries with type="outer"
   
  bpt<-rasterToPoints(regionb,fun=function(x){x>0},spatial=F)
  
  #select the most southerly point on the boundary
  distN<- rasterize(as.data.frame(bpt[order(bpt[,"y"]),1:2])[1,],regionb,field=2,update=T)

  #calculate stepwise distance around boundary from the south to all other points
  distN2<-gridDistance(x=distN, origin=2, omit=c(0,NA,NaN))
  
  distN<- distN2/cellStats(distN2, stat='max', na.rm=TRUE)
  
  bpt<-rasterToPoints(distN,spatial=F)
  
  #select points closest to halfway round in either direction
  ew<- bpt[order(abs(bpt[,"layer"]-0.5)),][1:16,]
  
  #of those, select the most westerly for new start
  distE<- rasterize(as.data.frame(ew[order(ew[,"x"]),1:2])[1,],regionb,field=2,update=T)
  
  #calculate stepwise distance around boundary from the west to all other points
  distE2<-gridDistance(x=distE, origin=2, omit=c(0,NA,NaN))
  
  distE<- distE2/cellStats(distE2, stat='max', na.rm=TRUE)
  
  ##now combine northerliness and easterliness in a data frame of points
  
  bpt<-as.data.frame(bpt)
  
  names(bpt)[3]<-"distN"
  bpt<-merge(bpt,as.data.frame(rasterToPoints(distE,spatial=F)))
  names(bpt)[4]<-"distE"
  
  #categorise the boundary based on relative distances N and E
  bpt$sector<- paste( c("S","N")[1+(bpt$distN>0.5)],c("W","E")[1+(bpt$distE>0.5)])
  
  bpt$sector[(bpt$distN>0.875)]<-"N" 
  bpt$sector[(bpt$distN<0.125)]<-"S"
  bpt$sector[(bpt$distE>0.875)]<-"E"
  bpt$sector[(bpt$distE<0.125)]<-"W"
  
  bpt$sector<-factor(bpt$sector)
  
  bpt
}


speednflow<-function(apt,sourcept,targetpt,R,disp,cellside,speedonly=TRUE){
  #condatis overall speed and the flow through each habitat cell, after raster data have been converted to points
  #apt should be a data frame with columns c("x","y","cover") - proportional coverage, no zeros included
  #sourcept and targetpt should be data frames with columns c("x","y")
  #cellside is the assumed width of each raster cell
  #disp is mean dispersal distance
  #use the same units for everything, either m or km, do not mix
  #R is number of emigrants produced per unit area (a sq-m if units m, a sq-km if units km), per generation
  
  len<-dim(apt)[1]
  dm<-dist(apt[,c("x","y")])
  alpha<-2/disp#alpha should be 2/dispersal distance
  
  norm<-R*alpha^2/2/pi*cellside^4
  
  Cfree<- norm*outer(apt$cover, apt$cover, "*")*exp(-alpha*as.matrix(dm))
  diag(Cfree)<-0
  
  Cin<- norm*apt$cover*
    rowSums(exp(-alpha*sqrt( outer(apt[,"x"],sourcept[,"x"],'-')^2 + outer(apt[,"y"],sourcept[,"y"],'-')^2 ) ))
  Cout<- norm*apt$cover*
    rowSums(exp(-alpha*sqrt( outer(apt[,"x"],targetpt[,"x"],'-')^2 + outer(apt[,"y"],targetpt[,"y"],'-')^2 ) ))
  
  
  M0 <- diag(Cin + Cout + rowSums(Cfree)) - Cfree
  w <- Cin - Cout
  ###voltage (note that I set voltage of target to -1 rather than 0. This seems to give more stable results)
  v0 <- solve(M0, w, tol = exp(-255))
  
  ###overall conductance
  I0 <- (v0 + 1) * Cout
  I1 <- (1 - v0) * Cin
  cond <- (sum(I0) + sum(I1))/4 #4 to make overall voltage difference 1 rather than 2
  
  if(speedonly){return(list(speed=cond))}else{
    
    ###flow by cell
    cur <- Cfree * outer(v0, v0, "-")
    flo <- apply(cur, 1, function(x) {
      sum(abs(x))
    })
    flo <- ((flo)/2 + I0 + I1)/2
    
    return(list(speed=cond,flow=flo))
  }#end if
}#end function




