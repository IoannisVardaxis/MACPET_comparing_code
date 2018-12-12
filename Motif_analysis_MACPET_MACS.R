########################################################################################################
######## script for running 5 runs of motivs for MACPET and MACS (and the different parameters for MACS)
########################################################################################################
# Important: One has to specify a save directory and a data directory, the current module only gives information about how the models were ran.
library(rGADEM)
library(MACPET)
library(BSgenome.Hsapiens.UCSC.hg19)#suggested in the MANUALL
library(MotIV)
library(plyr)
library(intervals)
########################################################################################################
######################################################################################################## MACPET
########################################################################################################
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZZ"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results/S3_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACPET_psfitData")
psfitData=MACPET_psfitData
data=metadata(psfitData)$Peaks.Info
data$point.estimate=round(data$Peak.Summit)
data$PE.s100=data$point.estimate-100#start 100
data$PE.e100=data$point.estimate+100# end 100
cat("\n")
data=data[with(data,order(FDR,decreasing=F)),]
data=data[1:5000,]#take the needed subset top 5000
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")
# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$Chrom)
    gadem.MACPET<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACPET,file="gadem.MACPET")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACPET)
    MotIV.MACPET <- motifMatch(motifs)#match motifs with the gadem
    save(MotIV.MACPET,file="MotIV.MACPET")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("ESR1"))
    Expected.motif.filtered=filter(MotIV.MACPET,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACPET)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(Chrom),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$Chrom))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }

                return(x.i)

            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }

        return(x)

    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACPET")
    cat(paste("Statistic MACPET A1: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}
rm(list=ls())
########################################################################################################
######################################################################################################## MACS
########################################################################################################
# For A1:
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZZ"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACS_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACS.peaks.out_peaks.update")
data=MACS.peaks.out_peaks.update
cat("\n")
data=data[with(data,order(p.values.FDR,decreasing=F)),]
data=data[1:5000,]
data$point.estimate=data$start+data$summit
data$PE.s100=data$point.estimate-100
data$PE.e100=data$point.estimate+100
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")
# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$chr)
    gadem.MACS<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACS,file="gadem.MACS")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACS)
    MotIV.MACS <- motifMatch(motifs)
    save(MotIV.MACS,file="MotIV.MACS")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("ESR1"))
    Expected.motif.filtered=filter(MotIV.MACS,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACS)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(chr),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$chr))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }

                return(x.i)

            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }

        return(x)

    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACS")
    cat(paste("Statistic MACS A1: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}
rm(list=ls())

########################################################################################################
######################################################################################################## MACPET
########################################################################################################
# For A3:
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAC"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results/S3_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACPET_psfitData")
psfitData=MACPET_psfitData
data=metadata(psfitData)$Peaks.Info
data$point.estimate=round(data$Peak.Summit)
data$PE.s100=data$point.estimate-100#start 100
data$PE.e100=data$point.estimate+100# end 100
data=data[with(data,order(FDR,decreasing=F)),]
data=data[1:5000,]#take the needed subset
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")
# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$Chrom)
    gadem.MACPET<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACPET,file="gadem.MACPET")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACPET)
    MotIV.MACPET <- motifMatch(motifs)#match motifs with the gadem
    save(MotIV.MACPET,file="MotIV.MACPET")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("CTCF"))
    Expected.motif.filtered=filter(MotIV.MACPET,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACPET)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(Chrom),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$Chrom))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }

                return(x.i)

            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }

        return(x)

    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACPET")
    cat(paste("Statistic MACPET A3: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}
rm(list=ls())

########################################################################################################
######################################################################################################## MACS
########################################################################################################
# For A3:
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAC"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACS_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACS.peaks.out_peaks.update")
data=MACS.peaks.out_peaks.update
# data=subset(data,p.values.FDR<1e-5)
# Sign=nrow(data)
# cat(Sign)
# cat("\n")
data=data[with(data,order(p.values.FDR,decreasing=F)),]
data=data[1:5000,]
data$point.estimate=data$start+data$summit
data$PE.s100=data$point.estimate-100
data$PE.e100=data$point.estimate+100
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")
# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$chr)
    gadem.MACS<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACS,file="gadem.MACS")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACS)
    MotIV.MACS <- motifMatch(motifs)
    save(MotIV.MACS,file="MotIV.MACS")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("CTCF"))
    Expected.motif.filtered=filter(MotIV.MACS,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACS)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(chr),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$chr))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }

                return(x.i)

            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }

        return(x)

    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACS")
    cat(paste("Statistic MACS A3: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}
rm(list=ls())

#######################################################################################################
####################################################################################################### MACPET
#######################################################################################################
# For A2:
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAD"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results/S3_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACPET_psfitData")
psfitData=MACPET_psfitData
data=metadata(psfitData)$Peaks.Info
data$point.estimate=round(data$Peak.Summit)
data$PE.s100=data$point.estimate-100#start 100
data$PE.e100=data$point.estimate+100# end 100
data=data[with(data,order(FDR,decreasing=F)),]
data=data[1:5000,]#take the needed subset
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")
# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$Chrom)
    gadem.MACPET<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACPET,file="gadem.MACPET")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACPET)
    MotIV.MACPET <- motifMatch(motifs)#match motifs with the gadem
    save(MotIV.MACPET,file="MotIV.MACPET")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("CTCF"))
    Expected.motif.filtered=filter(MotIV.MACPET,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACPET)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(Chrom),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$Chrom))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }

                return(x.i)

            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }

        return(x)

    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACPET")
    cat(paste("Statistic MACPET A2: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}

########################################################################################################
######################################################################################################## MACS
########################################################################################################
# For A2:
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAD"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACS_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACS.peaks.out_peaks.update")
data=MACS.peaks.out_peaks.update
data=subset(data,p.values.FDR<1e-5)
cat("\n")
data=data[with(data,order(p.values.FDR,decreasing=F)),]
data=data[1:5000,]
data$point.estimate=data$start+data$summit
data$PE.s100=data$point.estimate-100
data$PE.e100=data$point.estimate+100
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")

# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$chr)
    gadem.MACS<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACS,file="gadem.MACS")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACS)
    MotIV.MACS <- motifMatch(motifs)
    save(MotIV.MACS,file="MotIV.MACS")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("CTCF"))
    Expected.motif.filtered=filter(MotIV.MACS,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACS)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(chr),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$chr))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }

                return(x.i)

            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }

        return(x)

    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACS")
    cat(paste("Statistic MACS A2: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}
rm(list=ls())

########################################################################################################
########################################################################################################
########################################################################################################
#                From here you have MACS with the different parameters MACS2,MACS3,MACS4

########################################################################################################
######################################################################################################## MACS 2
########################################################################################################
#-------------------------------------------
#-------------------------------------------For A1 MACS2:
#-------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZZ"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACS2_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACS.peaks.out_peaks.update")
data=MACS.peaks.out_peaks.update
data=data[with(data,order(p.values.FDR,decreasing=F)),]
data=data[1:5000,]
data$point.estimate=data$start+data$summit
data$PE.s100=data$point.estimate-100
data$PE.e100=data$point.estimate+100
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")

# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$chr)
    gadem.MACS<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACS,file="gadem.MACS")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACS)
    MotIV.MACS <- motifMatch(motifs)
    save(MotIV.MACS,file="MotIV.MACS")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("ESR1"))
    Expected.motif.filtered=filter(MotIV.MACS,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACS)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(chr),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$chr))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }

                return(x.i)

            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }

        return(x)

    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACS")
    cat(paste("Statistic MACS2 Experiment_ENCSR000BZZ: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}
rm(list=ls())

#-------------------------------------------
#-------------------------------------------For A2 MACS2:
#-------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAD"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACS2_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACS.peaks.out_peaks.update")
data=MACS.peaks.out_peaks.update
data=data[with(data,order(p.values.FDR,decreasing=F)),]
data=data[1:5000,]
data$point.estimate=data$start+data$summit
data$PE.s100=data$point.estimate-100
data$PE.e100=data$point.estimate+100
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")

# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$chr)
    gadem.MACS<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACS,file="gadem.MACS")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACS)
    MotIV.MACS <- motifMatch(motifs)
    save(MotIV.MACS,file="MotIV.MACS")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("CTCF"))
    Expected.motif.filtered=filter(MotIV.MACS,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACS)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(chr),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$chr))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }
                
                return(x.i)
                
            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }
        
        return(x)
        
    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACS")
    cat(paste("Statistic MACS2 Experiment_ENCSR000CAD: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}
rm(list=ls())

#-------------------------------------------
#-------------------------------------------For A3 MACS2:
#-------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAC"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACS2_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACS.peaks.out_peaks.update")
data=MACS.peaks.out_peaks.update
data=data[with(data,order(p.values.FDR,decreasing=F)),]
data=data[1:5000,]
data$point.estimate=data$start+data$summit
data$PE.s100=data$point.estimate-100
data$PE.e100=data$point.estimate+100
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")

# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$chr)
    gadem.MACS<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACS,file="gadem.MACS")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACS)
    MotIV.MACS <- motifMatch(motifs)
    save(MotIV.MACS,file="MotIV.MACS")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("CTCF"))
    Expected.motif.filtered=filter(MotIV.MACS,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACS)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(chr),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$chr))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }
                
                return(x.i)
                
            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }
        
        return(x)
        
    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACS")
    cat(paste("Statistic MACS2 Experiment_ENCSR000CAC: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}
rm(list=ls())

########################################################################################################
######################################################################################################## MACS 3
########################################################################################################
#-------------------------------------------
#-------------------------------------------For A1 MACS3:
#-------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZZ"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACS3_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACS.peaks.out_peaks.update")
data=MACS.peaks.out_peaks.update
data=data[with(data,order(p.values.FDR,decreasing=F)),]
data=data[1:5000,]
data$point.estimate=data$start+data$summit
data$PE.s100=data$point.estimate-100
data$PE.e100=data$point.estimate+100
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")

# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$chr)
    gadem.MACS<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACS,file="gadem.MACS")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACS)
    MotIV.MACS <- motifMatch(motifs)
    save(MotIV.MACS,file="MotIV.MACS")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("ESR1"))
    Expected.motif.filtered=filter(MotIV.MACS,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACS)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(chr),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$chr))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }

                return(x.i)

            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }

        return(x)

    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACS")
    cat(paste("Statistic MACS3 Experiment_ENCSR000BZZ: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}
rm(list=ls())

#-------------------------------------------
#-------------------------------------------For A2 MACS3:
#-------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAD"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACS3_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACS.peaks.out_peaks.update")
data=MACS.peaks.out_peaks.update
data=data[with(data,order(p.values.FDR,decreasing=F)),]
data=data[1:5000,]
data$point.estimate=data$start+data$summit
data$PE.s100=data$point.estimate-100
data$PE.e100=data$point.estimate+100
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")

# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$chr)
    gadem.MACS<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACS,file="gadem.MACS")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACS)
    MotIV.MACS <- motifMatch(motifs)
    save(MotIV.MACS,file="MotIV.MACS")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("CTCF"))
    Expected.motif.filtered=filter(MotIV.MACS,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACS)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(chr),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$chr))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }

                return(x.i)

            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }

        return(x)

    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACS")
    cat(paste("Statistic MACS3 Experiment_ENCSR000CAD: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}
rm(list=ls())

#-------------------------------------------
#-------------------------------------------For A3 MACS3:
#-------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAC"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACS3_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACS.peaks.out_peaks.update")
data=MACS.peaks.out_peaks.update
data=data[with(data,order(p.values.FDR,decreasing=F)),]
data=data[1:5000,]
data$point.estimate=data$start+data$summit
data$PE.s100=data$point.estimate-100
data$PE.e100=data$point.estimate+100
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")

# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$chr)
    gadem.MACS<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACS,file="gadem.MACS")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACS)
    MotIV.MACS <- motifMatch(motifs)
    save(MotIV.MACS,file="MotIV.MACS")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("CTCF"))
    Expected.motif.filtered=filter(MotIV.MACS,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACS)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(chr),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$chr))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }

                return(x.i)

            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }

        return(x)

    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACS")
    cat(paste("Statistic MACS3 Experiment_ENCSR000CAC: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}
rm(list=ls())

########################################################################################################
######################################################################################################## MACS 4
########################################################################################################
#-------------------------------------------
#-------------------------------------------For A1 MACS4:
#-------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZZ"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACS4_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACS.peaks.out_peaks.update")
data=MACS.peaks.out_peaks.update
data=data[with(data,order(p.values.FDR,decreasing=F)),]
data=data[1:5000,]
data$point.estimate=data$start+data$summit
data$PE.s100=data$point.estimate-100
data$PE.e100=data$point.estimate+100
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")

# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$chr)
    gadem.MACS<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACS,file="gadem.MACS")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACS)
    MotIV.MACS <- motifMatch(motifs)
    save(MotIV.MACS,file="MotIV.MACS")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("ESR1"))
    Expected.motif.filtered=filter(MotIV.MACS,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACS)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(chr),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$chr))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }

                return(x.i)

            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }

        return(x)

    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACS")
    cat(paste("Statistic MACS4 Experiment_ENCSR000BZZ: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}
rm(list=ls())

#-------------------------------------------
#-------------------------------------------For A2 MACS4:
#-------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAD"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACS4_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACS.peaks.out_peaks.update")
data=MACS.peaks.out_peaks.update
data=data[with(data,order(p.values.FDR,decreasing=F)),]
data=data[1:5000,]
data$point.estimate=data$start+data$summit
data$PE.s100=data$point.estimate-100
data$PE.e100=data$point.estimate+100
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")

# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$chr)
    gadem.MACS<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACS,file="gadem.MACS")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACS)
    MotIV.MACS <- motifMatch(motifs)
    save(MotIV.MACS,file="MotIV.MACS")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("CTCF"))
    Expected.motif.filtered=filter(MotIV.MACS,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACS)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(chr),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$chr))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }

                return(x.i)

            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }

        return(x)

    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACS")
    cat(paste("Statistic MACS4 Experiment_ENCSR000CAD: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}
rm(list=ls())

#-------------------------------------------
#-------------------------------------------For A3 MACS4:
#-------------------------------------------
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAC"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACS4_results")
setwd(DIR)
path=system.file(package="MotIV")
jaspar=readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
load("MACS.peaks.out_peaks.update")
data=MACS.peaks.out_peaks.update
data=data[with(data,order(p.values.FDR,decreasing=F)),]
data=data[1:5000,]
data$point.estimate=data$start+data$summit
data$PE.s100=data$point.estimate-100
data$PE.e100=data$point.estimate+100
if(!dir.exists("Motifs")) dir.create("Motifs")
setwd("Motifs")
if(!dir.exists("Motifs.5000")) dir.create("Motifs.5000")
setwd("Motifs.5000")

# loop:
for(i in c(1:5)){
    if(!dir.exists(paste("RUN_",i,sep=""))) dir.create(paste("RUN_",i,sep=""))
    setwd(paste("RUN_",i,sep=""))
    #################################################################
    rgBED<-IRanges(start=data$PE.s100,end=data$PE.e100)
    Sequences<-RangedData(rgBED,space=data$chr)
    gadem.MACS<-GADEM(Sequences,verbose=1,genome=Hsapiens,seed=NULL)#note the genome has to be the same as the one from the data, here hg19
    save(gadem.MACS,file="gadem.MACS")
    jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
    motifs <- getPWM(gadem.MACS)
    MotIV.MACS <- motifMatch(motifs)
    save(MotIV.MACS,file="MotIV.MACS")
    dir.create("Expected")
    setwd("Expected")
    Expected.motif=setFilter(tfname=c("CTCF"))
    Expected.motif.filtered=filter(MotIV.MACS,Expected.motif,exact=F)
    Expected.motif.filtered.RD=exportAsRangedData(Expected.motif.filtered,gadem.MACS)
    Expected.motif.filtered.RD.df=as.data.frame(Expected.motif.filtered.RD)
    save(Expected.motif.filtered.RD.df,file="Expected.motif.filtered.RD.df")
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    Expected.motif.filtered.RD.df$Mid=round(rowMeans(Expected.motif.filtered.RD.df[,c("start","end")]))
    data.motif.overlap.expected=ddply(data,.(chr),function(x,Expected.motif.filtered.RD.df){
        #--------initialize:
        x$Expected.motif=0
        x$Dist.to.expected.motif=NA
        #get subset of the chromosome:
        Expected.motif.filtered.RD.df.chr=subset(Expected.motif.filtered.RD.df,space==unique(x$chr))
        if(nrow(Expected.motif.filtered.RD.df.chr)!=0){
            #----create intervals:
            PeakInt=intervals::Intervals(x[,c("PE.s100","PE.e100")],closed=c(T,T))
            MotifInt=intervals::Intervals(Expected.motif.filtered.RD.df.chr[,c("start","end")],closed=c(T,T))
            #----find overlaps:
            Ovlp=intervals::interval_included(PeakInt,MotifInt)
            #for each peak, keep the closest motif if any:
            x=ldply(1:length(Ovlp),function(i,x,Ovlp,Expected.motif.filtered.RD.df.chr){
                #take the Ovlp.i:
                Ovlp.i=Ovlp[[i]]
                #take x.i:
                x.i=x[i,]
                if(length(Ovlp.i)!=0){
                    x.i$Expected.motif=1#since it has a motif
                    #take Expected.motif.filtered.RD.df.chr.i:
                    Expected.motif.filtered.RD.df.chr.i=Expected.motif.filtered.RD.df.chr[Ovlp.i,]
                    Dist.to.expected.motif.i=abs(x.i$point.estimate-Expected.motif.filtered.RD.df.chr.i$Mid)
                    x.i$Dist.to.expected.motif=min(Dist.to.expected.motif.i)
                }

                return(x.i)

            },x=x,Expected.motif.filtered.RD.df.chr=Expected.motif.filtered.RD.df.chr,Ovlp=Ovlp)
        }

        return(x)

    },Expected.motif.filtered.RD.df=Expected.motif.filtered.RD.df,.progress="text")
    save(data.motif.overlap.expected,file="data.motif.overlap.expected.MACS")
    cat(paste("Statistic MACS4 Experiment_ENCSR000CAC: ","RUN ",i,"\n"))
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif))))#peaks with relevant motivs
    cat("\n")
    cat(nrow(subset(data.motif.overlap.expected,!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50)))#peaks with relevant motivs with dist<50
    cat("\n")
    #################################################################
    setwd("..")#go back
    setwd("..")#go back
}
rm(list=ls())

