##############################################################
##############################################################
## Script title: adding habitat lowish function
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
## Purpose of script: possible method of adding patches
##                    in "between-y" locations
##                    
## Warning: although much of this code has potential to be re-used 
## in projects other than the one for which this script was intended,
## it was not written with broad re-usability in mind. Only
## minimal comments are included for guidance 
##############################################################
##############################################################

# load libraries
library(raster)
library(spdep)
library(rgeos)
library(rgdal)



addhabitat.lowish<-  function(options){
  

  #extract the options
  regionfname<- options$regionfname # label to use for study region
  hablabel<- options$hablabel # label to use for study habitat
  habtext<- options$habtext #e.g. "semi.natural.grassland",#what the focal habitat is called in the sizes lookup table
  datapath<- options$datapath #file path for writing outputs, e.g. getwd()
  findfhabitatreg<- options$findfhabitatreg # file location for baseline habitat raster
  findhk2500<- options$findhk2500 # object name for longt-distance Hanski connectivity calculated on baseline habitat raster by function 'hanski'
  fhabitatbasics<- options$fhabitatbasics  # list of landscape stats as created by 'makehabitat' function
  addinglevels<- options$addinglevels # vector of proportions of habitat to add, e.g. c(0.01,0.05)
  resn<-options$resn # resolution of all rasters
  sizesfname<-  options$sizesfname #location of csv file of patch size distributions, e.g. "sizes.csv"#sizesfname should include any necessary path
  
  #load data about starting landscape

  hk2500<-get(findhk2500)
  fhabitatreg <- raster(findfhabitatreg)
  
  #load patch area distribution
  sizes<-read.csv(sizesfname)
  #create data frame to return, to check basic stats
  verify<-expand.grid(add=addinglevels,method=c("flattish"))
  
  #############zones of particularly low connectivity slope#############
  
  consl<-terrain(hk2500$conr, opt=c('slope'), unit='degrees',neighbors=8)
  
  #plot(consl)
  #plot(hk2500$conr)
  
  qq<- quantile(hk2500$conr,seq(0,1,0.03125))
  hkz<-cut(hk2500$conr,breaks=unique(qq))
  
  #plot(hkz)
  
  
  z1<-zonal(consl,hkz,fun=function(x,...){quantile(x,probs=0.1,na.rm=T)})
  slthresh<-subs(hkz,as.data.frame(z1))
  
  
  #mask out habitat too
  lowish<- mask((slthresh-consl)>0, mask=fhabitatreg, maskvalue=1, updatevalue=0)
  
  
  #########################################rest of function##########################
  #load patch area distribution
  #sizes<-read.csv(options$sizesfname)#factor-doesn't work
  
  sizes<-read.csv(as.character(options$sizesfname))
  
  
  ###set max sample size 
  
  propnonhab<- (1-fhabitatbasics$proparea)*(fhabitatbasics$regarea/
                                              (ncell(fhabitatreg)*resn*resn/1000/1000))
  
  propzone<- cellStats(lowish,sum)/ncell(lowish)
  
  
  ###create distributions of patch diameter and patch area
  didist<- rep(sqrt(sizes$npix),times= round(sizes[,habtext]*1000))
  ardist<- didist^2*resn*resn/1000/1000
  
  #hist(didist)
  #hist(ardist)
  
  #max sample determined as 10* that needed at the max level of existing habitat area, given the mean area of each patch and loss of pixels randomly outside non-habitat and 
  #  maxsamp<- round(10* fhabitatbasics$totarea/mean(ardist) * max(addinglevels) /propnonhab )
  
  #will select randomly within zone so don't need much more than target
  maxsamp<- round(2* fhabitatbasics$totarea/mean(ardist) * max(addinglevels) /propzone )
  
  
  chosen<- sampleRandom(fhabitatreg, size=maxsamp, na.rm=F,  
                        cells=TRUE, rowcol=FALSE, xy=T, sp=FALSE, asRaster=FALSE)
  chosen<-as.data.frame(chosen)
  
  #dim(chosen)
  
  nhchosen<-chosen[!is.na(chosen[,4]) & chosen[,4]==0,]
  
  
  ##remove picked locations outside zone
  
  nhchosen$zone<- extract(x=lowish,y=nhchosen[,c("x","y")], method = "simple")
  
  # summary(nhchosen)
  
  nhchosen<-nhchosen[!is.na(nhchosen[,"zone"]) & nhchosen[,"zone"]==1,]
  
  rm(chosen)
  
  #dim(nhchosen)
  
  ####assign areas to candidates
  
  nhchosen$di<- sample(didist, dim(nhchosen)[1], replace =T)
  #  base::
  nhchosen$patchid <- 2:(1+ dim(nhchosen)[1]) #the existing habitat will have id 1
  
  
  #summary(nhchosen)
  
  #####  #add cells to create bigger patches (around cell centres)
  
  morecells<-nhchosen
  centrecells<-nhchosen
  
  maxring<- (sqrt(max(sizes$npix))-1)/2
  
  for(i in c(1:maxring)){#this will be the no. of cells to count away in each direction
    
    side<-centrecells[(centrecells$di -1)/2 >= i,]
    
    if(dim(side)[1]>0){
      sideall<- side
      sideall$x<- sideall$x - (i*resn)
      sideall$y<- sideall$y - (i*resn)
      
      for(j in (-i+1):i){
        
        #left side
        side2<- side
        side2$x<- side2$x - (i*resn)
        side2$y<- side2$y + (j*resn)
        
        sideall<- rbind(sideall,side2)
        
        #right side
        side2<- side
        side2$y<- side2$y - (i*resn)
        side2$x<- side2$x + (j*resn)
        
        sideall<- rbind(sideall,side2)
        
        #right side
        side2<- side
        side2$x<- side2$x + (i*resn)
        side2$y<- side2$y + (j*resn)
        
        sideall<- rbind(sideall,side2)
        
        #top side
        side2<- side
        side2$y<- side2$y + (i*resn)
        side2$x<- side2$x + (j*resn)
        
        if (j<i){#this to prevent last, top-right cell being duplicated
          sideall<- rbind(sideall,side2)
        }#end if
      }#end j loop: filling 4 sides of square
      morecells<-rbind(morecells,sideall)
    }}#end if and i loop
  
  cat(paste0("\n", 'Cells added', "\n", Sys.time(), "\n"))
  
  #candidates not overlapping each other, favour lowest cell number for consistency
  
  morecells2<-aggregate(morecells[,c(1,4)],by=morecells[,c("x","y")],FUN=min)
  
  #summary(morecells2) 
  
  morecells2<-merge(morecells2,morecells,all.x=T,all.y=F)
  
  #  summary(morecells2)  
  
  
  #######candidates not overlapping existing habitat
  
  nonhab <- subs(x=fhabitatreg,y=data.frame(id=1, v=NA), subsWithNA=F)
  
  morecells<-rasterize(x=morecells2[,c("x","y")],y=nonhab, field=as.matrix(morecells2$patchid),update=TRUE,updateValue='!NA')
  
  #par(mfrow=c(1,1))
  #plot(fhabitatreg)
  #plot(morecells,col=rainbow(200)) 
  #plot(lowish,col=c(NA,"blue"),add=T) 
  
  #######turn back to points for selection
  
  morecells2<-rasterToPoints(morecells,fun=function(x){x>0})
  
  #######aggregate by id
  
  morecells2<-as.data.frame(morecells2)
  names(morecells2)[3]<-"patchid"
  
  nhchosen<-nhchosen[nhchosen$patchid %in% morecells2[,"patchid"],] #so,only centres of those patches that have some area left
  
  
  nhchosen$arealeft<- as.vector(table(morecells2$patchid)[as.character(nhchosen$patchid)])  *resn*resn/1000/1000
  
  #don't need connectivity values
  
  ####### place thresholds and get random subsets
  
  addingtargets<- fhabitatbasics$totarea* addinglevels
  maxartarget<- fhabitatbasics$totarea* max(addinglevels,na.rm=T)
  
  #sum(nhchosen$arealeft)
  #maxartarget# there is enough
  
  
  nhchosen<-nhchosen[order(nhchosen$patchid),] 
  
  
  ##add buffer to area target to prevent consistently undershooting
  addingtargets<- addingtargets+ mean(nhchosen$arealeft,na.rm=T)/2
  
  
  ###merge different thresholded versions with original habitat
  
  nhchosen$cumarea<-cumsum(nhchosen$arealeft)
  
  nhchosen$cumarea<-cut(nhchosen$cumarea,breaks=c(0,addingtargets,max(nhchosen$cumarea)))
  
  #summary(nhchosen)
  

  for(i in 1:length(addingtargets)){
    
    pickedids<- nhchosen$patchid[as.numeric(nhchosen$cumarea)<=i]
    
    morecells<-rasterize(x=morecells2[morecells2$patchid %in% pickedids ,c("x","y")],y=fhabitatreg, field=1,update=TRUE,updateValue='all')
    
    writeRaster(morecells, filename=paste0(datapath, "/", paste(regionfname,hablabel,"flattish.plus",addinglevels[i]*100,"tif",sep=".")),overwrite=TRUE)
    verify$nadded[verify$add==addinglevels[i] & verify$method=="flattish"] <- length(pickedids)
    
    verify$areaadded[verify$add==addinglevels[i] & verify$method=="flattish"] <- aadded
    
 
    
  }
  
verify
  
}#end the habitat adding funtion