###########################################################################################
################################ Create Figures for the article: ########################## This is very messy, but it only for creating the figures used in the article.
###########################################################################################
# Important: One has to specify a save directory and a data directory, the current module only gives information about how the models were ran.
library(MACPET)
library(ggplot2)
library(gridExtra)#side by side plots
library(ggpubr)
library(cowplot)
library(grid)
library(VennDiagram)
library(ggforce)
#################################################################################################################
################################ Figures for the cut-off: #################################
#################################################################################################################
scientific_10 <- function(x) {
    y=parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
    if(x[1]==0) y[1]=0
    return(y)
}#for scientific values
# load data for each dataset and make the figure:
#----------------
#----------------A1: ENCSR000BZZ
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZZ"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
DIRsaveA1=file.path(SaveDir,"ENCSR000BZZ")
if(!dir.exists(DIRsaveA1)) dir.create(DIRsaveA1)
load(file.path(DIR,"S2_results/MACPET_pselfData"))
pselfData=MACPET_pselfData
load(file.path(DIR,"S2_results/MACPET_pintraData"))
pintraData=MACPET_pintraData
DataA1=c(pselfData,pintraData)
DataA1$Dist=InteractionSet::pairdist(DataA1,type="span")
MaxSpanA1=max(DataA1$Dist)
MinSpanA1=min(DataA1$Dist)
CutOfA1=metadata(pselfData)$MaxSize
SeqSpanA1=seq(from=0,to=MaxSpanA1,by=100)
SeqSpanA1[which(SeqSpanA1==max(SeqSpanA1))]=MaxSpanA1
IntSpanA1=cbind(SeqSpanA1[-length(SeqSpanA1)],SeqSpanA1[-1])
IntSpanA1=intervals::Intervals(IntSpanA1,closed=rep(TRUE,2))
DataIntA1=intervals::Intervals(cbind(DataA1$Dist,DataA1$Dist),closed=rep(TRUE,2))
IntInclA1=intervals::interval_included(IntSpanA1,DataIntA1)
# save("IntInclA1",file=file.path(SaveDir,"IntInclA1"))
# load(file.path(SaveDir,"IntInclA1"))
SpanDFA1=data.frame(Size=SeqSpanA1[-1],Freq=lengths(IntInclA1))
SpanDFA1=subset(SpanDFA1,Freq!=0)
SpanDFA1$logSize=log(SpanDFA1$Size)
save("SpanDFA1",file=file.path(DIRsaveA1,"SpanDFA1"))
save("CutOfA1",file=file.path(DIRsaveA1,"CutOfA1"))
# plot:
A1SelfCutPlot=
    ggplot2::ggplot(subset(SpanDFA1,Freq>3),ggplot2::aes(x=Size,y=Freq))+
    # add histogram:
    ggplot2::geom_rect(ggplot2::aes(xmin=subset(SpanDFA1,Freq>3)$Size-49,
                                    xmax=subset(SpanDFA1,Freq>3)$Size+49,
                                    ymin=0,ymax=subset(SpanDFA1,Freq>3)$Freq),
                       size=0.3,fill="grey69",color="black")+
    # add line:
    ggplot2::geom_line(color="blue",size=0.3)+
    # add cut-off line:
    ggplot2::geom_vline(xintercept=CutOfA1,
                        linetype="dashed",color="red",size=0.3)+
    # change scale x and remove empty space from the left:
    ggplot2::scale_x_log10(labels=function(x) round(log(x)), expand = c(0, 0),
                           breaks=exp(pretty(log(subset(SpanDFA1,Freq>3)$Size),n=10)))+
    # add titles
    ggplot2::ylab("")+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(label=scientific_10,expand = c(0,0))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    # add self cut:
    annotate("text", x=11*CutOfA1, y=sort(SpanDFA1$Freq,decreasing=T)[3], label= paste("Self-ligated cut-off at ",as.character(CutOfA1)," bp",sep=""),size=2.5)

A1SelfCutPlot
# save the plot as R object, then it can be loaded and create the figure
save("A1SelfCutPlot",file=file.path(DIRsaveA1,"A1SelfCutPlot"))

#----------------
#----------------A2: ENCSR000CAD
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAD"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
DIRsaveA2=file.path(SaveDir,"ENCSR000CAD")
if(!dir.exists(DIRsaveA2)) dir.create(DIRsaveA2)
load(file.path(DIR,"S2_results/MACPET_pselfData"))
pselfData=MACPET_pselfData
load(file.path(DIR,"S2_results/MACPET_pintraData"))
pintraData=MACPET_pintraData
DataA2=c(pselfData,pintraData)
DataA2$Dist=InteractionSet::pairdist(DataA2,type="span")
MaxSpanA2=max(DataA2$Dist)
MinSpanA2=min(DataA2$Dist)
CutOfA2=metadata(pselfData)$MaxSize
SeqSpanA2=seq(from=0,to=MaxSpanA2,by=100)
SeqSpanA2[which(SeqSpanA2==max(SeqSpanA2))]=MaxSpanA2
IntSpanA2=cbind(SeqSpanA2[-length(SeqSpanA2)],SeqSpanA2[-1])
IntSpanA2=intervals::Intervals(IntSpanA2,closed=rep(TRUE,2))
DataIntA2=intervals::Intervals(cbind(DataA2$Dist,DataA2$Dist),closed=rep(TRUE,2))
IntInclA2=intervals::interval_included(IntSpanA2,DataIntA2)
# save("IntInclA2",file=file.path(SaveDir,"IntInclA2"))
# load(file.path(SaveDir,"IntInclA2"))
SpanDFA2=data.frame(Size=SeqSpanA2[-1],Freq=lengths(IntInclA2))
SpanDFA2=subset(SpanDFA2,Freq!=0)
SpanDFA2$logSize=log(SpanDFA2$Size)
save("SpanDFA2",file=file.path(DIRsaveA2,"SpanDFA2"))
save("CutOfA2",file=file.path(DIRsaveA2,"CutOfA2"))
# plot:
A2SelfCutPlot=
    ggplot2::ggplot(subset(SpanDFA2,Freq>3),ggplot2::aes(x=Size,y=Freq))+
    # add histogram:
    ggplot2::geom_rect(ggplot2::aes(xmin=subset(SpanDFA2,Freq>3)$Size-49,
                                    xmax=subset(SpanDFA2,Freq>3)$Size+49,
                                    ymin=0,ymax=subset(SpanDFA2,Freq>3)$Freq),
                       size=0.3,fill="grey69",color="black")+
    # add line:
    ggplot2::geom_line(color="blue",size=0.3)+
    # add cut-off line:
    ggplot2::geom_vline(xintercept=CutOfA2,
                        linetype="dashed",color="red",size=0.3)+
    # change scale x and remove empty space from the left:
    ggplot2::scale_x_log10(labels=function(x) round(log(x)), expand = c(0, 0),
                           breaks=exp(pretty(log(subset(SpanDFA2,Freq>3)$Size),n=10)))+
    # add titles
    ggplot2::ylab("Frequency")+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(label=scientific_10,expand = c(0, 0))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    # add self cut:
    annotate("text", x=11*CutOfA2, y=sort(SpanDFA2$Freq,decreasing=T)[2], label= paste("Self-ligated cut-off at ",as.character(CutOfA2)," bp",sep=""),size=2.5)


A2SelfCutPlot
save("A2SelfCutPlot",file=file.path(DIRsaveA2,"A2SelfCutPlot"))

#----------------
#----------------A3: ENCSR000CAC
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAC"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
DIRsaveA3=file.path(SaveDir,"ENCSR000CAC")
if(!dir.exists(DIRsaveA3)) dir.create(DIRsaveA3)
load(file.path(DIR,"S2_results/MACPET_pselfData"))
pselfData=MACPET_pselfData
load(file.path(DIR,"S2_results/MACPET_pintraData"))
pintraData=MACPET_pintraData
DataA3=c(pselfData,pintraData)
DataA3$Dist=InteractionSet::pairdist(DataA3,type="span")
MaxSpanA3=max(DataA3$Dist)
MinSpanA3=min(DataA3$Dist)
CutOfA3=metadata(pselfData)$MaxSize
SeqSpanA3=seq(from=0,to=MaxSpanA3,by=100)
SeqSpanA3[which(SeqSpanA3==max(SeqSpanA3))]=MaxSpanA3
IntSpanA3=cbind(SeqSpanA3[-length(SeqSpanA3)],SeqSpanA3[-1])
IntSpanA3=intervals::Intervals(IntSpanA3,closed=rep(TRUE,2))
DataIntA3=intervals::Intervals(cbind(DataA3$Dist,DataA3$Dist),closed=rep(TRUE,2))
IntInclA3=intervals::interval_included(IntSpanA3,DataIntA3)
# save("IntInclA3",file=file.path(SaveDir,"IntInclA3"))
# load(file.path(SaveDir,"IntInclA3"))
SpanDFA3=data.frame(Size=SeqSpanA3[-1],Freq=lengths(IntInclA3))
SpanDFA3=subset(SpanDFA3,Freq!=0)
SpanDFA3$logSize=log(SpanDFA3$Size)
save("SpanDFA3",file=file.path(DIRsaveA3,"SpanDFA3"))
save("CutOfA3",file=file.path(DIRsaveA3,"CutOfA3"))
# plot:
A3SelfCutPlot=
    ggplot2::ggplot(subset(SpanDFA3,Freq>3),ggplot2::aes(x=Size,y=Freq))+
    # add histogram:
    ggplot2::geom_rect(ggplot2::aes(xmin=subset(SpanDFA3,Freq>3)$Size-49,
                                    xmax=subset(SpanDFA3,Freq>3)$Size+49,
                                    ymin=0,ymax=subset(SpanDFA3,Freq>3)$Freq),
                       size=0.3,fill="grey69",color="black")+
    # add line:
    ggplot2::geom_line(color="blue",size=0.3)+
    # add cut-off line:
    ggplot2::geom_vline(xintercept=CutOfA3,
                        linetype="dashed",color="red",size=0.3)+
    # change scale x and remove empty space from the left:
    ggplot2::scale_x_log10(labels=function(x) round(log(x)), expand = c(0, 0),
                           breaks=exp(pretty(log(subset(SpanDFA3,Freq>3)$Size),n=10)))+
    # add titles
    ggplot2::xlab("Sorted PET sizes in log10-scale")+
    ggplot2::ylab("")+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(label=scientific_10,expand = c(0, 0))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    # add self cut:
    annotate("text", x=10*CutOfA3, y=sort(SpanDFA3$Freq,decreasing=T)[2], label= paste("Self-ligated cut-off at ",as.character(CutOfA3)," bp",sep=""),size=2.5)

A3SelfCutPlot
save("A3SelfCutPlot",file=file.path(DIRsaveA3,"A3SelfCutPlot"))


#----------------
#----------------A4: ENCSR000FDD
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000FDD"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
DIRsaveA4=file.path(SaveDir,"ENCSR000FDD")
if(!dir.exists(DIRsaveA4)) dir.create(DIRsaveA4)
load(file.path(DIR,"S2_results/MACPET_pselfData"))
pselfData=MACPET_pselfData
load(file.path(DIR,"S2_results/MACPET_pintraData"))
pintraData=MACPET_pintraData
DataA4=c(pselfData,pintraData)
DataA4$Dist=InteractionSet::pairdist(DataA4,type="span")
MaxSpanA4=max(DataA4$Dist)
MinSpanA4=min(DataA4$Dist)
CutOfA4=metadata(pselfData)$MaxSize
SeqSpanA4=seq(from=0,to=MaxSpanA4,by=100)
SeqSpanA4[which(SeqSpanA4==max(SeqSpanA4))]=MaxSpanA4
IntSpanA4=cbind(SeqSpanA4[-length(SeqSpanA4)],SeqSpanA4[-1])
IntSpanA4=intervals::Intervals(IntSpanA4,closed=rep(TRUE,2))
DataIntA4=intervals::Intervals(cbind(DataA4$Dist,DataA4$Dist),closed=rep(TRUE,2))
IntInclA4=intervals::interval_included(IntSpanA4,DataIntA4)
# save("IntInclA4",file=file.path(SaveDir,"IntInclA4"))
# load(file.path(SaveDir,"IntInclA4"))
SpanDFA4=data.frame(Size=SeqSpanA4[-1],Freq=lengths(IntInclA4))
SpanDFA4=subset(SpanDFA4,Freq!=0)
SpanDFA4$logSize=log(SpanDFA4$Size)
save("SpanDFA4",file=file.path(DIRsaveA4,"SpanDFA4"))
save("CutOfA4",file=file.path(DIRsaveA4,"CutOfA4"))
# plot:
A4SelfCutPlot=
    ggplot2::ggplot(subset(SpanDFA4,Freq>3),ggplot2::aes(x=Size,y=Freq))+
    # add histogram:
    ggplot2::geom_rect(ggplot2::aes(xmin=subset(SpanDFA4,Freq>3)$Size-49,
                                    xmax=subset(SpanDFA4,Freq>3)$Size+49,
                                    ymin=0,ymax=subset(SpanDFA4,Freq>3)$Freq),
                       size=0.3,fill="grey69",color="black")+
    # add line:
    ggplot2::geom_line(color="blue",size=0.3)+
    # add cut-off line:
    ggplot2::geom_vline(xintercept=CutOfA4,
                        linetype="dashed",color="red",size=0.3)+
    # change scale x and remove empty space from the left:
    ggplot2::scale_x_log10(labels=function(x) round(log(x)), expand = c(0, 0),
                           breaks=exp(pretty(log(subset(SpanDFA4,Freq>3)$Size),n=10)))+
    # add titles
    ggplot2::ylab("")+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(label=scientific_10,expand = c(0,0))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    # add self cut:
    annotate("text", x=11*CutOfA4, y=sort(SpanDFA4$Freq,decreasing=T)[2], label= paste("Self-ligated cut-off at ",as.character(CutOfA4)," bp",sep=""),size=2.5)

A4SelfCutPlot
save("A4SelfCutPlot",file=file.path(DIRsaveA4,"A4SelfCutPlot"))


#----------------
#----------------A5: ENCSR000FDG
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000FDG"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
DIRsaveA5=file.path(SaveDir,"ENCSR000FDG")
if(!dir.exists(DIRsaveA5)) dir.create(DIRsaveA5)
load(file.path(DIR,"S2_results/MACPET_pselfData"))
pselfData=MACPET_pselfData
load(file.path(DIR,"S2_results/MACPET_pintraData"))
pintraData=MACPET_pintraData
DataA5=c(pselfData,pintraData)
DataA5$Dist=InteractionSet::pairdist(DataA5,type="span")
MaxSpanA5=max(DataA5$Dist)
MinSpanA5=min(DataA5$Dist)
CutOfA5=metadata(pselfData)$MaxSize
SeqSpanA5=seq(from=0,to=MaxSpanA5,by=100)
SeqSpanA5[which(SeqSpanA5==max(SeqSpanA5))]=MaxSpanA5
IntSpanA5=cbind(SeqSpanA5[-length(SeqSpanA5)],SeqSpanA5[-1])
IntSpanA5=intervals::Intervals(IntSpanA5,closed=rep(TRUE,2))
DataIntA5=intervals::Intervals(cbind(DataA5$Dist,DataA5$Dist),closed=rep(TRUE,2))
IntInclA5=intervals::interval_included(IntSpanA5,DataIntA5)
# save("IntInclA4",file=file.path(SaveDir,"IntInclA4"))
# load(file.path(SaveDir,"IntInclA4"))
SpanDFA5=data.frame(Size=SeqSpanA5[-1],Freq=lengths(IntInclA5))
SpanDFA5=subset(SpanDFA5,Freq!=0)
SpanDFA5$logSize=log(SpanDFA5$Size)
save("SpanDFA5",file=file.path(DIRsaveA5,"SpanDFA5"))
save("CutOfA5",file=file.path(DIRsaveA5,"CutOfA5"))
# plot:
A5SelfCutPlot=
    ggplot2::ggplot(subset(SpanDFA5,Freq>3),ggplot2::aes(x=Size,y=Freq))+
    # add histogram:
    ggplot2::geom_rect(ggplot2::aes(xmin=subset(SpanDFA5,Freq>3)$Size-49,
                                    xmax=subset(SpanDFA5,Freq>3)$Size+49,
                                    ymin=0,ymax=subset(SpanDFA5,Freq>3)$Freq),
                       size=0.3,fill="grey69",color="black")+
    # add line:
    ggplot2::geom_line(color="blue",size=0.3)+
    # add cut-off line:
    ggplot2::geom_vline(xintercept=CutOfA5,
                        linetype="dashed",color="red",size=0.3)+
    # change scale x and remove empty space from the left:
    ggplot2::scale_x_log10(labels=function(x) round(log(x)), expand = c(0, 0),
                           breaks=exp(pretty(log(subset(SpanDFA5,Freq>3)$Size),n=10)))+
    # add titles
    ggplot2::ylab("Frequency")+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(label=scientific_10,expand = c(0, 0))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    # add self cut:
    annotate("text", x=11*CutOfA5, y=sort(SpanDFA5$Freq,decreasing=T)[2], label= paste("Self-ligated cut-off at ",as.character(CutOfA5)," bp",sep=""),size=2.5)

A5SelfCutPlot
save("A5SelfCutPlot",file=file.path(DIRsaveA5,"A5SelfCutPlot"))

#----------------
#----------------A6: ENCSR000BZY
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZY"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
DIRsaveA6=file.path(SaveDir,"ENCSR000BZY")
if(!dir.exists(DIRsaveA6)) dir.create(DIRsaveA6)
load(file.path(DIR,"S2_results/MACPET_pselfData"))
pselfData=MACPET_pselfData
load(file.path(DIR,"S2_results/MACPET_pintraData"))
pintraData=MACPET_pintraData
DataA6=c(pselfData,pintraData)
DataA6$Dist=InteractionSet::pairdist(DataA6,type="span")
MaxSpanA6=max(DataA6$Dist)
MinSpanA6=min(DataA6$Dist)
CutOfA6=metadata(pselfData)$MaxSize
SeqSpanA6=seq(from=0,to=MaxSpanA6,by=100)
SeqSpanA6[which(SeqSpanA6==max(SeqSpanA6))]=MaxSpanA6
IntSpanA6=cbind(SeqSpanA6[-length(SeqSpanA6)],SeqSpanA6[-1])
IntSpanA6=intervals::Intervals(IntSpanA6,closed=rep(TRUE,2))
DataIntA6=intervals::Intervals(cbind(DataA6$Dist,DataA6$Dist),closed=rep(TRUE,2))
IntInclA6=intervals::interval_included(IntSpanA6,DataIntA6)
SpanDFA6=data.frame(Size=SeqSpanA6[-1],Freq=lengths(IntInclA6))
SpanDFA6=subset(SpanDFA6,Freq!=0)
SpanDFA6$logSize=log(SpanDFA6$Size)
save("SpanDFA6",file=file.path(DIRsaveA6,"SpanDFA6"))
save("CutOfA6",file=file.path(DIRsaveA6,"CutOfA6"))
# plot:
A6SelfCutPlot=
    ggplot2::ggplot(subset(SpanDFA6,Freq>3),ggplot2::aes(x=Size,y=Freq))+
    # add histogram:
    ggplot2::geom_rect(ggplot2::aes(xmin=subset(SpanDFA6,Freq>3)$Size-49,
                                    xmax=subset(SpanDFA6,Freq>3)$Size+49,
                                    ymin=0,ymax=subset(SpanDFA6,Freq>3)$Freq),
                       size=0.3,fill="grey69",color="black")+
    # add line:
    ggplot2::geom_line(color="blue",size=0.3)+
    # add cut-off line:
    ggplot2::geom_vline(xintercept=CutOfA6,
                        linetype="dashed",color="red",size=0.3)+
    # change scale x and remove empty space from the left:
    ggplot2::scale_x_log10(labels=function(x) round(log(x)), expand = c(0, 0),
                           breaks=exp(pretty(log(subset(SpanDFA6,Freq>3)$Size),n=10)))+
    # add titles
    ggplot2::xlab("Sorted PET sizes in log10-scale")+
    ggplot2::ylab("")+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(label=scientific_10,expand = c(0, 0))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    # add self cut:
    annotate("text", x=11*CutOfA6, y=sort(SpanDFA6$Freq,decreasing=T)[2], label= paste("Self-ligated cut-off at ",as.character(CutOfA6)," bp",sep=""),size=2.5)

A6SelfCutPlot
save("A6SelfCutPlot",file=file.path(DIRsaveA6,"A6SelfCutPlot"))


#################################################################################################################
########################################### Figures for region: ###################################### Figure 2
#################################################################################################################
# The saving directory
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DirSaveAll=file.path(SaveDir,"Lever")
if(!dir.exists(DirSaveAll)) dir.create(DirSaveAll)

# load psfit data:
load("YOURPATH/ChIA_PET_Data/Experiment_ENCSR000CAD/REP1/MACPET_results/S3_results/MACPET_psfitData")
x=MACPET_psfitData
# Peaks.Info=S4Vectors::metadata(x)$Peaks.Info
# Peaks.Info=subset(Peaks.Info,FDR<1e-5)
# choose region:
Classification.Info=S4Vectors::metadata(x)$Classification.Info
xdf=data.frame(x)
Chrom=xdf$seqnames1
Chrom=as.character(Chrom)
Classification.Info$Chrom=Chrom[Classification.Info$MainIndex]
Classification.Info$RegCount=paste(Classification.Info$Region,"-",
                                   Classification.Info$Chrom,sep="")
MaxReg=table(Classification.Info$RegCount)
MaxReg=sort(MaxReg,decreasing=TRUE)
i=14#6,7
which.id.region=names(MaxReg)[i]
#take subset which will be plotted:
Classification.Info_i=subset(Classification.Info,
                           RegCount==which.id.region)
xsub=xdf[Classification.Info_i$MainIndex,]
xsub$Peak.ID=0
xsub$Peak.ID=Classification.Info_i$Peak.ID
xsub=xsub[,c("start1","end1","start2","end2","Peak.ID")]
#sort:
Tosort=which(xsub$start1>xsub$start2)
if(length(Tosort)>0){
    Start2new=xsub$start1[Tosort]
    End2new=xsub$end1[Tosort]
    strand2new=xsub$strand1[Tosort]
    xsub$start1[Tosort]=xsub$start2[Tosort]
    xsub$end1[Tosort]=xsub$end2[Tosort]
    xsub$strand1[Tosort]=xsub$strand2[Tosort]
    xsub$start2[Tosort]=Start2new
    xsub$end2[Tosort]=End2new
    xsub$strand2[Tosort]=strand2new
}
xsub$Dist=xsub$end2-xsub$start1+1
#take the  peak summits of the region:
Peaks.Info=S4Vectors::metadata(x)$Peaks.Info
Peaks.Info$RegCount=paste(Peaks.Info$Region,"-",
                          Peaks.Info$Chrom,sep="")
Peaks.Info=subset(Peaks.Info,RegCount==which.id.region)
Peaks.Info=subset(Peaks.Info,FDR<5e-2)
xsub$Peak.ID[which(!xsub$Peak.ID%in%Peaks.Info$Peak)]=0
# plot:
#plot region PETs, no strand info.
Dens=data.frame(Y=(xsub$start1+xsub$end2)/2,ymin=xsub$start1,
                ymax=xsub$end2,
                X=xsub$Dist,PeakID=xsub$Peak.ID)
# plot:
RegionPlot=ggplot2::ggplot(Dens, ggplot2::aes(x=X,y=Y,ymin=ymin,ymax=ymax,
                                       color=factor(PeakID)))+
    # scale of y and removing space beneath:
    ggplot2::scale_x_continuous(expand = c(0, 0),breaks=seq(from=0,to=max(Dens$X),by=200))+
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    # error bars for the PEts
    ggplot2::geom_errorbar(width = 35,show.legend = F,size=0.1)+
    ggplot2::coord_flip()+
    ggplot2::geom_hline(yintercept=Peaks.Info$Peak.Summit,
                        color="black",linetype="dashed",size=0.3)+
    ggplot2::ylab("Midpoints of PETs on the genome")+
    ggplot2::xlab("PET sizes")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=7,family="Times"),
          axis.title=element_text(size=9,family="Times"))




RegionPlot

# save:
ggsave(filename=file.path(DirSaveAll,"Figure_8.pdf"),
       plot=RegionPlot,
       dpi = 300,units="mm",width=170,height=112.5,device="pdf",family="Times")

#isws na kaneis kai ta densities, alla thes ta estimations


#################################################################################################################
################################ Figures for Peak overlaps and sizes: ###############################
#################################################################################################################
# The saving directory
#----------------
#----------------A1: ENCSR000BZZ
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZZ"
Repetition="REP1"
DIRsaveA1=file.path(SaveDir,"ENCSR000BZZ")
############ Load data MACPET ###############
MACPETpeaks=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results/MACPET_psfitData",sep="")# peaks data
load(MACPETpeaks)
MACPETpeaks=metadata(MACPET_psfitData)$Peaks.Info
nrow(subset(MACPETpeaks,FDR<0.05))
nrow(subset(MACPETpeaks,FDR<1e-5))
#make GRanges:
MACPETpeaksSign=subset(MACPETpeaks,FDR<0.05)
MACPETpeaksGRangesA1=GenomicRanges::GRanges(seqnames=MACPETpeaksSign$Chrom,
                                          ranges=IRanges::IRanges(start=MACPETpeaksSign$CIQ.Up.start,
                                                                  end=MACPETpeaksSign$CIQ.Down.end))
############ Load data MACS ###############
load(paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACS_results/MACS.peaks.out_peaks.update",sep=""))# peaks data)
nrow(subset(MACS.peaks.out_peaks.update,p.values.FDR<0.05))
nrow(subset(MACS.peaks.out_peaks.update,p.values.FDR<1e-5))
#make GRanges:
MACSpeaksSign=subset(MACS.peaks.out_peaks.update,p.values.FDR<0.05)
MACSpeaksGRangesA1=GenomicRanges::GRanges(seqnames=MACSpeaksSign$chr,
                                        ranges=IRanges::IRanges(start=MACSpeaksSign$start,
                                                                end=MACSpeaksSign$end))
############ Find overlaps ###############
# make calculations for creating the venn data:

TotPeaksMACPETA1=length(MACPETpeaksGRangesA1)
TotPeaksMACSA1=length(MACSpeaksGRangesA1)
TotPeaksA1=TotPeaksMACPETA1+TotPeaksMACSA1
TotCommonPeaksMACStoMACPETA1=InteractionSet::countOverlaps(MACSpeaksGRangesA1,MACPETpeaksGRangesA1)
TotPeaks_overlapA1=length(which(TotCommonPeaksMACStoMACPETA1!=0))
# change to 100 scale
TotPeaksMACPETA1_100=TotPeaksMACPETA1*100/TotPeaksA1
TotPeaksMACSA1_100=TotPeaksMACSA1*100/TotPeaksA1
TotPeaks_overlapA1_100=TotPeaks_overlapA1*100/TotPeaksA1

VennPeaksA1_data=data.frame(x0=c(0+TotPeaksMACPETA1_100/2+TotPeaks_overlapA1_100/2,#macpet
                                 100-TotPeaksMACSA1_100/2-TotPeaks_overlapA1_100/2,NA),#macs
                            y0=c(0.5,0.5,NA),
                            r=c(TotPeaksMACPETA1_100/2,TotPeaksMACSA1_100/2,NA),
                            Totals=c(TotPeaksMACPETA1,TotPeaksMACSA1,TotPeaks_overlapA1),
                            color=c("red","blue",NA),which=c("MACPET","MACS",NA))

save("VennPeaksA1_data",file=file.path(DIRsaveA1,"VennPeaksA1_data"))


VennPeaksA1_plot=ggplot()+geom_circle(data=VennPeaksA1_data[2,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="blue",alpha=0.6,size=0.2)+
    geom_circle(data=VennPeaksA1_data[1,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="red",alpha=0.6,size=0.2)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),axis.ticks=element_blank(),
          axis.text=element_blank(),axis.title=element_blank())+
    # add numbers for overlaps:
    annotate("text", x=15, y=0.5, label= as.character(VennPeaksA1_data$Totals[3]),size=2.5)+#common
    annotate("text", x=100-VennPeaksA1_data$r[2]-10, y=0.5, label= as.character(VennPeaksA1_data$Totals[2]-VennPeaksA1_data$Totals[3]),size=2.5)+#macs
    annotate("text", x=0+VennPeaksA1_data$r[1], y=22, label= as.character(VennPeaksA1_data$Totals[1]-VennPeaksA1_data$Totals[3]),size=2.5)+#MACPET
    # add line for the text:
    geom_line(data=data.frame(x=rep(0+VennPeaksA1_data$r[1],2),y=c(0,20)),aes(x,y),size=0.3)+
    # add groups:
    annotate("text", x=15, y=40, label= "MACPET",size=3)+#MACPET
    annotate("text", x=85, y=40, label= "MACS",size=3)#MACPET


VennPeaksA1_plot
save("VennPeaksA1_plot",file=file.path(DIRsaveA1,"VennPeaksA1_plot"))

#plot the strengths and sizes now:
PeakStrengthSize_A1_data=data.frame(StrengthTags=c(2*MACPETpeaksSign$Pets,MACSpeaksSign$tags),
                                Which=c(rep("MACPET",nrow(MACPETpeaksSign)),rep("MACS",nrow(MACSpeaksSign))),
                                Size=c(MACPETpeaksSign$CIQ.Peak.size,MACSpeaksSign$length))
save("PeakStrengthSize_A1_data",file=file.path(DIRsaveA1,"PeakStrengthSize_A1_data"))

PeakStrength_A1_plot=ggplot(PeakStrengthSize_A1_data,aes(StrengthTags,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6,size=0.5)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0.1),trans="log10",
                                breaks=round(exp(pretty(log(PeakStrengthSize_A1_data$StrengthTags),n=5))))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    xlab("")+ylab("Density")+
    # ylab("Frequency")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

PeakStrength_A1_plot
save("PeakStrength_A1_plot",file=file.path(DIRsaveA1,"PeakStrength_A1_plot"))

PeakSizes_A1_plot=ggplot(PeakStrengthSize_A1_data,aes(Size,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6,size=0.5)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0.1),
                                breaks=round((pretty((PeakStrengthSize_A1_data$Size),n=5))))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    ylab("Density")+
    xlab("")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

PeakSizes_A1_plot
save("PeakSizes_A1_plot",file=file.path(DIRsaveA1,"PeakSizes_A1_plot"))

#----------------
#----------------A2: ENCSR000CAD
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAD"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
DIRsaveA2=file.path(SaveDir,"ENCSR000CAD")
############ Load data MACPET ###############
MACPETpeaks=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results/MACPET_psfitData",sep="")# peaks data
load(MACPETpeaks)
MACPETpeaks=metadata(MACPET_psfitData)$Peaks.Info
nrow(subset(MACPETpeaks,FDR<0.05))
nrow(subset(MACPETpeaks,FDR<1e-5))
#make GRanges:
MACPETpeaksSign=subset(MACPETpeaks,FDR<0.05)
MACPETpeaksGRangesA2=GenomicRanges::GRanges(seqnames=MACPETpeaksSign$Chrom,
                                          ranges=IRanges::IRanges(start=MACPETpeaksSign$CIQ.Up.start,
                                                                  end=MACPETpeaksSign$CIQ.Down.end))
############ Load data MACS ###############
load(paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACS_results/MACS.peaks.out_peaks.update",sep=""))# peaks data)
nrow(subset(MACS.peaks.out_peaks.update,p.values.FDR<0.05))
nrow(subset(MACS.peaks.out_peaks.update,p.values.FDR<1e-5))
#make GRanges:
MACSpeaksSign=subset(MACS.peaks.out_peaks.update,p.values.FDR<0.05)
MACSpeaksGRangesA2=GenomicRanges::GRanges(seqnames=MACSpeaksSign$chr,
                                        ranges=IRanges::IRanges(start=MACSpeaksSign$start,
                                                                end=MACSpeaksSign$end))
############ Find overlaps ###############
# make calculations for creating the venn data:
TotPeaksMACPETA2=length(MACPETpeaksGRangesA2)
TotPeaksMACSA2=length(MACSpeaksGRangesA2)
TotPeaksA2=TotPeaksMACPETA2+TotPeaksMACSA2
TotCommonPeaksMACStoMACPETA2=InteractionSet::countOverlaps(MACSpeaksGRangesA2,MACPETpeaksGRangesA2)
TotPeaks_overlapA2=length(which(TotCommonPeaksMACStoMACPETA2!=0))
# change to 100 scale
TotPeaksMACPETA2_100=TotPeaksMACPETA2*100/TotPeaksA2
TotPeaksMACSA2_100=TotPeaksMACSA2*100/TotPeaksA2
TotPeaks_overlapA2_100=TotPeaks_overlapA2*100/TotPeaksA2

VennPeaksA2_data=data.frame(x0=c(0+TotPeaksMACPETA2_100/2+TotPeaks_overlapA2_100/2,#macpet
                                 100-TotPeaksMACSA2_100/2-TotPeaks_overlapA2_100/2,NA),#macs
                            y0=c(0.5,0.5,NA),
                            r=c(TotPeaksMACPETA2_100/2,TotPeaksMACSA2_100/2,NA),
                            Totals=c(TotPeaksMACPETA2,TotPeaksMACSA2,TotPeaks_overlapA2),
                            color=c("red","blue",NA),which=c("MACPET","MACS",NA))

save("VennPeaksA2_data",file=file.path(DIRsaveA2,"VennPeaksA2_data"))


VennPeaksA2_plot=ggplot()+geom_circle(data=VennPeaksA2_data[2,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="blue",alpha=0.6,size=0.2)+
    geom_circle(data=VennPeaksA2_data[1,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="red",alpha=0.6,size=0.2)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),axis.ticks=element_blank(),
          axis.text=element_blank(),axis.title=element_blank())+
    # add numbers for overlaps:
    annotate("text", x=25, y=0.5, label= as.character(VennPeaksA2_data$Totals[3]),size=2.5)+#common
    annotate("text", x=100-VennPeaksA2_data$r[2]-10, y=0.5, label= as.character(VennPeaksA2_data$Totals[2]-VennPeaksA2_data$Totals[3]),size=2.5)+#macs
    annotate("text", x=0+VennPeaksA2_data$r[1], y=22, label= as.character(VennPeaksA2_data$Totals[1]-VennPeaksA2_data$Totals[3]),size=2.5)+#MACPET
    # add line for the text:
    geom_line(data=data.frame(x=rep(0+VennPeaksA2_data$r[1],2),y=c(0,20)),aes(x,y),size=0.3)+
    # add groups:
    annotate("text", x=20, y=35, label= "MACPET",size=3)+#MACPET
    annotate("text", x=80, y=35, label= "MACS",size=3)#MACPET


VennPeaksA2_plot
save("VennPeaksA2_plot",file=file.path(DIRsaveA2,"VennPeaksA2_plot"))
# plot strength and size
PeakStrengthSize_A2_data=data.frame(StrengthTags=c(2*MACPETpeaksSign$Pets,MACSpeaksSign$tags),
                                    Which=c(rep("MACPET",nrow(MACPETpeaksSign)),rep("MACS",nrow(MACSpeaksSign))),
                                    Size=c(MACPETpeaksSign$CIQ.Peak.size,MACSpeaksSign$length))
save("PeakStrengthSize_A2_data",file=file.path(DIRsaveA2,"PeakStrengthSize_A2_data"))

PeakStrength_A2_plot=ggplot(PeakStrengthSize_A2_data,aes(StrengthTags,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6,size=0.5)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0.1),trans="log10",
                                breaks=round(exp(pretty(log(PeakStrengthSize_A2_data$StrengthTags),n=5))))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    xlab("Total Tags in log scale")+
    ylab("")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

PeakStrength_A2_plot
save("PeakStrength_A2_plot",file=file.path(DIRsaveA2,"PeakStrength_A2_plot"))

PeakSizes_A2_plot=ggplot(PeakStrengthSize_A2_data,aes(Size,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6,size=0.5)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0.1),
                                breaks=round((pretty((PeakStrengthSize_A2_data$Size),n=6))))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    xlab("Binding Site Size")+
    ylab("")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

PeakSizes_A2_plot
save("PeakSizes_A2_plot",file=file.path(DIRsaveA2,"PeakSizes_A2_plot"))

#----------------
#----------------A3: ENCSR000CAC
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAC"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
DIRsaveA3=file.path(SaveDir,"ENCSR000CAC")
############ Load data MACPET ###############
MACPETpeaks=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results/MACPET_psfitData",sep="")# peaks data
load(MACPETpeaks)
MACPETpeaks=metadata(MACPET_psfitData)$Peaks.Info
nrow(subset(MACPETpeaks,FDR<0.05))
nrow(subset(MACPETpeaks,FDR<1e-5))
#make GRanges:
MACPETpeaksSign=subset(MACPETpeaks,FDR<0.05)
MACPETpeaksGRangesA3=GenomicRanges::GRanges(seqnames=MACPETpeaksSign$Chrom,
                                            ranges=IRanges::IRanges(start=MACPETpeaksSign$CIQ.Up.start,
                                                                    end=MACPETpeaksSign$CIQ.Down.end))
############ Load data MACS ###############
load(paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACS_results/MACS.peaks.out_peaks.update",sep=""))# peaks data)
nrow(subset(MACS.peaks.out_peaks.update,p.values.FDR<1e-5))
nrow(subset(MACS.peaks.out_peaks.update,p.values.FDR<1e-5))
#make GRanges:
MACSpeaksSign=subset(MACS.peaks.out_peaks.update,p.values.FDR<0.05)
MACSpeaksGRangesA3=GenomicRanges::GRanges(seqnames=MACSpeaksSign$chr,
                                          ranges=IRanges::IRanges(start=MACSpeaksSign$start,
                                                                  end=MACSpeaksSign$end))
############ Find overlaps ###############
# make calculations for creating the venn data:
TotPeaksMACPETA3=length(MACPETpeaksGRangesA3)
TotPeaksMACSA3=length(MACSpeaksGRangesA3)
TotPeaksA3=TotPeaksMACPETA3+TotPeaksMACSA3
TotCommonPeaksMACStoMACPETA3=InteractionSet::countOverlaps(MACPETpeaksGRangesA3,MACSpeaksGRangesA3)
TotPeaks_overlapA3=length(which(TotCommonPeaksMACStoMACPETA3!=0))
# change to 100 scale
TotPeaksMACPETA3_100=TotPeaksMACPETA3*100/TotPeaksA3
TotPeaksMACSA3_100=TotPeaksMACSA3*100/TotPeaksA3
TotPeaks_overlapA3_100=TotPeaks_overlapA3*100/TotPeaksA3

VennPeaksA3_data=data.frame(x0=c(0+TotPeaksMACPETA3_100/2+TotPeaks_overlapA3_100/2,#macpet
                                 100-TotPeaksMACSA3_100/2-TotPeaks_overlapA3_100/2,NA),#macs
                            y0=c(0.5,0.5,NA),
                            r=c(TotPeaksMACPETA3_100/2,TotPeaksMACSA3_100/2,NA),
                            Totals=c(TotPeaksMACPETA3,TotPeaksMACSA3,TotPeaks_overlapA3),
                            color=c("red","blue",NA),which=c("MACPET","MACS",NA))

save("VennPeaksA3_data",file=file.path(DIRsaveA3,"VennPeaksA3_data"))


VennPeaksA3_plot=ggplot()+geom_circle(data=VennPeaksA3_data[2,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="blue",alpha=0.6,size=0.2)+
    geom_circle(data=VennPeaksA3_data[1,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="red",alpha=0.6,size=0.2)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),axis.ticks=element_blank(),
          axis.text=element_blank(),axis.title=element_blank())+
    # add numbers for overlaps:
    annotate("text", x=20, y=0.5, label= as.character(VennPeaksA3_data$Totals[3]),size=2.5)+#common
    annotate("text", x=100-VennPeaksA3_data$r[2]-10, y=0.5, label= as.character(VennPeaksA3_data$Totals[2]-VennPeaksA3_data$Totals[3]),size=2.5)+#macs
    annotate("text", x=0+VennPeaksA3_data$r[1], y=22, label= as.character(VennPeaksA3_data$Totals[1]-VennPeaksA3_data$Totals[3]),size=2.5)+#MACPET
    # add line for the text:
    geom_line(data=data.frame(x=rep(0+VennPeaksA3_data$r[1],2),y=c(0,20)),aes(x,y),size=0.3)+
    # add groups:
    annotate("text", x=17, y=37, label= "MACPET",size=3)+#MACPET
    annotate("text", x=83, y=37, label= "MACS",size=3)#MACPET


VennPeaksA3_plot
save("VennPeaksA3_plot",file=file.path(DIRsaveA3,"VennPeaksA3_plot"))

#plot the strengths and sizes now:
PeakStrengthSize_A3_data=data.frame(StrengthTags=c(2*MACPETpeaksSign$Pets,MACSpeaksSign$tags),
                                    Which=c(rep("MACPET",nrow(MACPETpeaksSign)),rep("MACS",nrow(MACSpeaksSign))),
                                    Size=c(MACPETpeaksSign$CIQ.Peak.size,MACSpeaksSign$length))
save("PeakStrengthSize_A3_data",file=file.path(DIRsaveA3,"PeakStrengthSize_A3_data"))

PeakStrength_A3_plot=ggplot(PeakStrengthSize_A3_data,aes(StrengthTags,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6,size=0.5)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0.1),trans="log10",
                                breaks=round(exp(pretty(log(PeakStrengthSize_A3_data$StrengthTags),n=5))))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    xlab("  ")+
    ylab("  ")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

PeakStrength_A3_plot
save("PeakStrength_A3_plot",file=file.path(DIRsaveA3,"PeakStrength_A3_plot"))

PeakSizes_A3_plot=ggplot(PeakStrengthSize_A3_data,aes(Size,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6,size=0.5)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0.1),
                                breaks=round((pretty((PeakStrengthSize_A3_data$Size),n=4))))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    xlab("  ")+
    ylab("  ")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

PeakSizes_A3_plot
save("PeakSizes_A3_plot",file=file.path(DIRsaveA3,"PeakSizes_A3_plot"))

#----------------
#----------------A4: ENCSR000FDD
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000FDD"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
DIRsaveA4=file.path(SaveDir,"ENCSR000FDD")
############ Load data MACPET ###############
MACPETpeaks=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results/MACPET_psfitData",sep="")# peaks data
load(MACPETpeaks)
MACPETpeaks=metadata(MACPET_psfitData)$Peaks.Info
nrow(subset(MACPETpeaks,FDR<0.05))
nrow(subset(MACPETpeaks,FDR<1e-5))
#make GRanges:
MACPETpeaksSign=subset(MACPETpeaks,FDR<0.05)
MACPETpeaksGRangesA4=GenomicRanges::GRanges(seqnames=MACPETpeaksSign$Chrom,
                                            ranges=IRanges::IRanges(start=MACPETpeaksSign$CIQ.Up.start,
                                                                    end=MACPETpeaksSign$CIQ.Down.end))
############ Load data MACS ###############
load(paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACS_results/MACS.peaks.out_peaks.update",sep=""))# peaks data)
nrow(subset(MACS.peaks.out_peaks.update,p.values.FDR<0.05))
nrow(subset(MACS.peaks.out_peaks.update,p.values.FDR<1e-5))
#make GRanges:
MACSpeaksSign=subset(MACS.peaks.out_peaks.update,p.values.FDR<0.05)
MACSpeaksGRangesA4=GenomicRanges::GRanges(seqnames=MACSpeaksSign$chr,
                                          ranges=IRanges::IRanges(start=MACSpeaksSign$start,
                                                                  end=MACSpeaksSign$end))
############ Find overlaps ###############
# make calculations for creating the venn data:
TotPeaksMACPETA4=length(MACPETpeaksGRangesA4)
TotPeaksMACSA4=length(MACSpeaksGRangesA4)
TotPeaksA4=TotPeaksMACPETA4+TotPeaksMACSA4
TotCommonPeaksMACStoMACPETA4=InteractionSet::countOverlaps(MACPETpeaksGRangesA4,MACSpeaksGRangesA4)
TotPeaks_overlapA4=length(which(TotCommonPeaksMACStoMACPETA4!=0))
# change to 100 scale
TotPeaksMACPETA4_100=TotPeaksMACPETA4*100/TotPeaksA4
TotPeaksMACSA4_100=TotPeaksMACSA4*100/TotPeaksA4
TotPeaks_overlapA4_100=TotPeaks_overlapA4*100/TotPeaksA4

VennPeaksA4_data=data.frame(x0=c(0+TotPeaksMACPETA4_100/2+TotPeaks_overlapA4_100/2,#macpet
                                 100-TotPeaksMACSA4_100/2-TotPeaks_overlapA4_100/2,NA),#macs
                            y0=c(0.5,0.5,NA),
                            r=c(TotPeaksMACPETA4_100/2,TotPeaksMACSA4_100/2,NA),
                            Totals=c(TotPeaksMACPETA4,TotPeaksMACSA4,TotPeaks_overlapA4),
                            color=c("red","blue",NA),which=c("MACPET","MACS",NA))

save("VennPeaksA4_data",file=file.path(DIRsaveA4,"VennPeaksA4_data"))


VennPeaksA4_plot=ggplot()+geom_circle(data=VennPeaksA4_data[2,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="blue",alpha=0.6,size=0.2)+
    geom_circle(data=VennPeaksA4_data[1,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="red",alpha=0.6,size=0.2)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),axis.ticks=element_blank(),
          axis.text=element_blank(),axis.title=element_blank())+
    # add numbers for overlaps:
    annotate("text", x=45, y=0.5, label= as.character(VennPeaksA4_data$Totals[3]),size=2.5)+#common
    annotate("text", x=100-VennPeaksA4_data$r[2], y=0.5, label= as.character(VennPeaksA4_data$Totals[2]-VennPeaksA4_data$Totals[3]),size=2.5)+#macs
    annotate("text", x=0+VennPeaksA4_data$r[1], y=22, label= as.character(VennPeaksA4_data$Totals[1]-VennPeaksA4_data$Totals[3]),size=2.5)+#MACPET
    # add line for the text:
    geom_line(data=data.frame(x=rep(0+VennPeaksA4_data$r[1],2),y=c(0,20)),aes(x,y),size=0.3)+
    # add groups:
    annotate("text", x=25, y=40, label= "MACPET",size=3)+#MACPET
    annotate("text", x=75, y=40, label= "MACS",size=3)#MACPET


VennPeaksA4_plot
save("VennPeaksA4_plot",file=file.path(DIRsaveA4,"VennPeaksA4_plot"))

# plot strength and size
PeakStrengthSize_A4_data=data.frame(StrengthTags=c(2*MACPETpeaksSign$Pets,MACSpeaksSign$tags),
                                    Which=c(rep("MACPET",nrow(MACPETpeaksSign)),rep("MACS",nrow(MACSpeaksSign))),
                                    Size=c(MACPETpeaksSign$CIQ.Peak.size,MACSpeaksSign$length))
save("PeakStrengthSize_A4_data",file=file.path(DIRsaveA4,"PeakStrengthSize_A4_data"))

PeakStrength_A4_plot=ggplot(PeakStrengthSize_A4_data,aes(StrengthTags,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6,size=0.5)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0.1),trans="log10",
                                breaks=round(exp(pretty(log(PeakStrengthSize_A4_data$StrengthTags),n=5))))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    xlab("")+
    ylab("Density")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

PeakStrength_A4_plot
save("PeakStrength_A4_plot",file=file.path(DIRsaveA4,"PeakStrength_A4_plot"))

PeakSizes_A4_plot=ggplot(PeakStrengthSize_A4_data,aes(Size,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6,size=0.5)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0.1),
                                breaks=round((pretty((PeakStrengthSize_A4_data$Size),n=3))))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    xlab("")+
    ylab("Density")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

PeakSizes_A4_plot
save("PeakSizes_A4_plot",file=file.path(DIRsaveA4,"PeakSizes_A4_plot"))

#----------------
#----------------A5: ENCSR000FDG or Not_specified_1
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000FDG"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
DIRsaveA5=file.path(SaveDir,"ENCSR000FDG")
############ Load data MACPET ###############
MACPETpeaks=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results/MACPET_psfitData",sep="")# peaks data
load(MACPETpeaks)
MACPETpeaks=metadata(MACPET_psfitData)$Peaks.Info
nrow(subset(MACPETpeaks,FDR<0.05))
nrow(subset(MACPETpeaks,FDR<1e-5))
#make GRanges:
MACPETpeaksSign=subset(MACPETpeaks,FDR<0.05)
MACPETpeaksGRangesA5=GenomicRanges::GRanges(seqnames=MACPETpeaksSign$Chrom,
                                            ranges=IRanges::IRanges(start=MACPETpeaksSign$CIQ.Up.start,
                                                                    end=MACPETpeaksSign$CIQ.Down.end))
############ Load data MACS ###############
load(paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACS_results/MACS.peaks.out_peaks.update",sep=""))# peaks data)
nrow(subset(MACS.peaks.out_peaks.update,p.values.FDR<0.05))
nrow(subset(MACS.peaks.out_peaks.update,p.values.FDR<1e-5))
#make GRanges:
MACSpeaksSign=subset(MACS.peaks.out_peaks.update,p.values.FDR<0.05)
MACSpeaksGRangesA5=GenomicRanges::GRanges(seqnames=MACSpeaksSign$chr,
                                          ranges=IRanges::IRanges(start=MACSpeaksSign$start,
                                                                  end=MACSpeaksSign$end))
############ Find overlaps ###############
# make calculations for creating the venn data:
TotPeaksMACPETA5=length(MACPETpeaksGRangesA5)
TotPeaksMACSA5=length(MACSpeaksGRangesA5)
TotPeaksA5=TotPeaksMACPETA5+TotPeaksMACSA5
TotCommonPeaksMACStoMACPETA5=InteractionSet::countOverlaps(MACPETpeaksGRangesA5,MACSpeaksGRangesA5)
TotPeaks_overlapA5=length(which(TotCommonPeaksMACStoMACPETA5!=0))
# change to 100 scale
TotPeaksMACPETA5_100=TotPeaksMACPETA5*100/TotPeaksA5
TotPeaksMACSA5_100=TotPeaksMACSA5*100/TotPeaksA5
TotPeaks_overlapA5_100=TotPeaks_overlapA5*100/TotPeaksA5

VennPeaksA5_data=data.frame(x0=c(0+TotPeaksMACPETA5_100/2+TotPeaks_overlapA5_100/2,#macpet
                                 100-TotPeaksMACSA5_100/2-TotPeaks_overlapA5_100/2,NA),#macs
                            y0=c(0.5,0.5,NA),
                            r=c(TotPeaksMACPETA5_100/2,TotPeaksMACSA5_100/2,NA),
                            Totals=c(TotPeaksMACPETA5,TotPeaksMACSA5,TotPeaks_overlapA5),
                            color=c("red","blue",NA),which=c("MACPET","MACS",NA))

save("VennPeaksA5_data",file=file.path(DIRsaveA5,"VennPeaksA5_data"))


VennPeaksA5_plot=ggplot()+geom_circle(data=VennPeaksA5_data[2,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="blue",alpha=0.6,size=0.2)+
    geom_circle(data=VennPeaksA5_data[1,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="red",alpha=0.6,size=0.2)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),axis.ticks=element_blank(),
          axis.text=element_blank(),axis.title=element_blank())+
    # add numbers for overlaps:
    annotate("text", x=32, y=0.5, label= as.character(VennPeaksA5_data$Totals[3]),size=2.5)+#common
    annotate("text", x=100-VennPeaksA5_data$r[2]-10, y=0.5, label= as.character(VennPeaksA5_data$Totals[2]-VennPeaksA5_data$Totals[3]),size=2.5)+#macs
    annotate("text", x=0+VennPeaksA5_data$r[1], y=22, label= as.character(VennPeaksA5_data$Totals[1]-VennPeaksA5_data$Totals[3]),size=2.5)+#MACPET
    # add line for the text:
    geom_line(data=data.frame(x=rep(0+VennPeaksA5_data$r[1],2),y=c(0,20)),aes(x,y),size=0.3)+
    # add groups:
    annotate("text", x=25, y=40, label= "MACPET",size=3)+#MACPET
    annotate("text", x=75, y=40, label= "MACS",size=3)#MACPET


VennPeaksA5_plot
save("VennPeaksA5_plot",file=file.path(DIRsaveA5,"VennPeaksA5_plot"))

# plot strength and size
PeakStrengthSize_A5_data=data.frame(StrengthTags=c(2*MACPETpeaksSign$Pets,MACSpeaksSign$tags),
                                    Which=c(rep("MACPET",nrow(MACPETpeaksSign)),rep("MACS",nrow(MACSpeaksSign))),
                                    Size=c(MACPETpeaksSign$CIQ.Peak.size,MACSpeaksSign$length))
save("PeakStrengthSize_A5_data",file=file.path(DIRsaveA5,"PeakStrengthSize_A5_data"))

PeakStrength_A5_plot=ggplot(PeakStrengthSize_A5_data,aes(StrengthTags,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6,size=0.5)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0.1),trans="log10",
                                breaks=round(exp(pretty(log(PeakStrengthSize_A5_data$StrengthTags),n=5))))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    xlab("Total Tags in log scale")+
    ylab("")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

PeakStrength_A5_plot
save("PeakStrength_A5_plot",file=file.path(DIRsaveA5,"PeakStrength_A5_plot"))

PeakSizes_A5_plot=ggplot(PeakStrengthSize_A5_data,aes(Size,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6,size=0.5)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0.1),
                                breaks=round((pretty((PeakStrengthSize_A5_data$Size),n=5))))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    xlab("Binding Site Size")+
    ylab("")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

PeakSizes_A5_plot
save("PeakSizes_A5_plot",file=file.path(DIRsaveA5,"PeakSizes_A5_plot"))

#----------------
#----------------A6: ENCSR000BZY
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZY"
Repetition="REP1"
DIR=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results")
DIRsaveA6=file.path(SaveDir,"ENCSR000BZY")
############ Load data MACPET ###############
MACPETpeaks=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results/MACPET_psfitData",sep="")# peaks data
load(MACPETpeaks)
MACPETpeaks=metadata(MACPET_psfitData)$Peaks.Info
nrow(subset(MACPETpeaks,FDR<0.05))
nrow(subset(MACPETpeaks,FDR<1e-5))
#make GRanges:
MACPETpeaksSign=subset(MACPETpeaks,FDR<0.05)
MACPETpeaksGRangesA6=GenomicRanges::GRanges(seqnames=MACPETpeaksSign$Chrom,
                                            ranges=IRanges::IRanges(start=MACPETpeaksSign$CIQ.Up.start,
                                                                    end=MACPETpeaksSign$CIQ.Down.end))
############ Load data MACS ###############
load(paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACS_results/MACS.peaks.out_peaks.update",sep=""))# peaks data)
nrow(subset(MACS.peaks.out_peaks.update,p.values.FDR<0.05))
nrow(subset(MACS.peaks.out_peaks.update,p.values.FDR<1e-5))
#make GRanges:
MACSpeaksSign=subset(MACS.peaks.out_peaks.update,p.values.FDR<0.05)
MACSpeaksGRangesA6=GenomicRanges::GRanges(seqnames=MACSpeaksSign$chr,
                                          ranges=IRanges::IRanges(start=MACSpeaksSign$start,
                                                                  end=MACSpeaksSign$end))
############ Find overlaps ###############
# make calculations for creating the venn data:
TotPeaksMACPETA6=length(MACPETpeaksGRangesA6)
TotPeaksMACSA6=length(MACSpeaksGRangesA6)
TotPeaksA6=TotPeaksMACPETA6+TotPeaksMACSA6
TotCommonPeaksMACStoMACPETA6=InteractionSet::countOverlaps(MACSpeaksGRangesA6,MACPETpeaksGRangesA6)
TotPeaks_overlapA6=length(which(TotCommonPeaksMACStoMACPETA6!=0))
# change to 100 scale
TotPeaksMACPETA6_100=TotPeaksMACPETA6*100/TotPeaksA6
TotPeaksMACSA6_100=TotPeaksMACSA6*100/TotPeaksA6
TotPeaks_overlapA6_100=TotPeaks_overlapA6*100/TotPeaksA6

VennPeaksA6_data=data.frame(x0=c(0+TotPeaksMACPETA6_100/2+TotPeaks_overlapA6_100/2,#macpet
                                 100-TotPeaksMACSA6_100/2-TotPeaks_overlapA6_100/2,NA),#macs
                            y0=c(0.5,0.5,NA),
                            r=c(TotPeaksMACPETA6_100/2,TotPeaksMACSA6_100/2,NA),
                            Totals=c(TotPeaksMACPETA6,TotPeaksMACSA6,TotPeaks_overlapA6),
                            color=c("red","blue",NA),which=c("MACPET","MACS",NA))

save("VennPeaksA6_data",file=file.path(DIRsaveA6,"VennPeaksA6_data"))


VennPeaksA6_plot=ggplot()+geom_circle(data=VennPeaksA6_data[2,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="blue",alpha=0.6,size=0.2)+
    geom_circle(data=VennPeaksA6_data[1,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="red",alpha=0.6,size=0.2)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),axis.ticks=element_blank(),
          axis.text=element_blank(),axis.title=element_blank())+
    # add numbers for overlaps:
    annotate("text", x=18, y=0.5, label= as.character(VennPeaksA6_data$Totals[3]),size=2.5)+#common
    annotate("text", x=100-VennPeaksA6_data$r[2]-10, y=0.5, label= as.character(VennPeaksA6_data$Totals[2]-VennPeaksA6_data$Totals[3]),size=2.5)+#macs
    annotate("text", x=0+VennPeaksA6_data$r[1], y=22, label= as.character(VennPeaksA6_data$Totals[1]-VennPeaksA6_data$Totals[3]),size=2.5)+#MACPET
    # add line for the text:
    geom_line(data=data.frame(x=rep(0+VennPeaksA6_data$r[1],2),y=c(0,20)),aes(x,y),size=0.3)+
    # add groups:
    annotate("text", x=15, y=37, label= "MACPET",size=3)+#MACPET
    annotate("text", x=85, y=37, label= "MACS",size=3)#MACPET


VennPeaksA6_plot
save("VennPeaksA6_plot",file=file.path(DIRsaveA6,"VennPeaksA6_plot"))

# plot strength and size
PeakStrengthSize_A6_data=data.frame(StrengthTags=c(2*MACPETpeaksSign$Pets,MACSpeaksSign$tags),
                                    Which=c(rep("MACPET",nrow(MACPETpeaksSign)),rep("MACS",nrow(MACSpeaksSign))),
                                    Size=c(MACPETpeaksSign$CIQ.Peak.size,MACSpeaksSign$length))
save("PeakStrengthSize_A6_data",file=file.path(DIRsaveA6,"PeakStrengthSize_A6_data"))

PeakStrength_A6_plot=ggplot(PeakStrengthSize_A6_data,aes(StrengthTags,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6,size=0.5)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0.1),trans="log10",
                                breaks=round(exp(pretty(log(PeakStrengthSize_A6_data$StrengthTags),n=4))))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    xlab("")+
    ylab("")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

PeakStrength_A6_plot
save("PeakStrength_A6_plot",file=file.path(DIRsaveA6,"PeakStrength_A6_plot"))

PeakSizes_A6_plot=ggplot(PeakStrengthSize_A6_data,aes(Size,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6,size=0.5)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0.1),
                                breaks=round((pretty((PeakStrengthSize_A6_data$Size),n=3))))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    xlab("")+
    ylab("")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

PeakSizes_A6_plot
save("PeakSizes_A6_plot",file=file.path(DIRsaveA6,"PeakSizes_A6_plot"))



#################################################################################################################
################################ Figures for MO, SE, FDR: ######################################
#################################################################################################################
# The saving directory
#----------------
#----------------A1: ENCSR000BZZ
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
############ Load data MACPET ###############
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZZ"
Repetition="REP1"
DIRsaveA1=file.path(SaveDir,"ENCSR000BZZ")
#save in the external drive:
DIR=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results",sep="")# the data directory
setwd(DIR)
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACPET.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACPET")
    data.MACPET.run_i=data.motif.overlap.expected
    data.MACPET.run_i$Run=i
    data.MACPET.Allruns=rbind(data.MACPET.Allruns,data.MACPET.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACPET.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACPET.run_1=subset(data.MACPET.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACPET.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACPET.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACPET.run_2=subset(data.MACPET.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACPET.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACPET.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACPET.run_3=subset(data.MACPET.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACPET.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACPET.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACPET.run_4=subset(data.MACPET.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACPET.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACPET.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACPET.run_5=subset(data.MACPET.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACPET.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACPET.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACPET.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                      Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                      Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                      Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                      Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                      Spatial.error.run_1=Spatial.error.run_1,
                                                      Spatial.error.run_2=Spatial.error.run_2,
                                                      Spatial.error.run_3=Spatial.error.run_3,
                                                      Spatial.error.run_4=Spatial.error.run_4,
                                                      Spatial.error.run_5=Spatial.error.run_5,
                                                      Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                      Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                      Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                  tot.motifs.run_2/tot.bs,
                                                                                  tot.motifs.run_3/tot.bs,
                                                                                  tot.motifs.run_4/tot.bs,
                                                                                  tot.motifs.run_5/tot.bs)),
                                                      Spatial.error.min=min(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                      Spatial.error.max=max(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                      Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                                Spatial.error.run_2,
                                                                                Spatial.error.run_3,
                                                                                Spatial.error.run_4,
                                                                                Spatial.error.run_5)),
                                                      tot.bs=x,Which="MACPET")

    # merge:
    data.MACPET.occurrence.spatial.error=rbind(data.MACPET.occurrence.spatial.error,
                                               data.MACPET.occurrence.spatial.error.x)

}
setwd("../../../..")
getwd()
############ Load data MACS ###############
setwd("MACS_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS.run_i=data.motif.overlap.expected
    data.MACS.run_i$Run=i
    data.MACS.Allruns=rbind(data.MACS.Allruns,data.MACS.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                    Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                    Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                    Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                    Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                    Spatial.error.run_1=Spatial.error.run_1,
                                                    Spatial.error.run_2=Spatial.error.run_2,
                                                    Spatial.error.run_3=Spatial.error.run_3,
                                                    Spatial.error.run_4=Spatial.error.run_4,
                                                    Spatial.error.run_5=Spatial.error.run_5,
                                                    Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                              tot.motifs.run_2/tot.bs,
                                                                              tot.motifs.run_3/tot.bs,
                                                                              tot.motifs.run_4/tot.bs,
                                                                              tot.motifs.run_5/tot.bs)),
                                                    Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                              tot.motifs.run_2/tot.bs,
                                                                              tot.motifs.run_3/tot.bs,
                                                                              tot.motifs.run_4/tot.bs,
                                                                              tot.motifs.run_5/tot.bs)),
                                                    Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                    Spatial.error.min=min(c(Spatial.error.run_1,
                                                                            Spatial.error.run_2,
                                                                            Spatial.error.run_3,
                                                                            Spatial.error.run_4,
                                                                            Spatial.error.run_5)),
                                                    Spatial.error.max=max(c(Spatial.error.run_1,
                                                                            Spatial.error.run_2,
                                                                            Spatial.error.run_3,
                                                                            Spatial.error.run_4,
                                                                            Spatial.error.run_5)),
                                                    Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                    tot.bs=x,Which="MACS")

    # merge:
    data.MACS.occurrence.spatial.error=rbind(data.MACS.occurrence.spatial.error,
                                             data.MACS.occurrence.spatial.error.x)

}
setwd("../../..")
getwd()
# setwd("MotifComparizons")
############ Create the whole data ###############
Motif.occurence.allA1=rbind(data.MACPET.occurrence.spatial.error,data.MACS.occurrence.spatial.error)
save("Motif.occurence.allA1",file=file.path(DIRsaveA1,"Motif.occurence.allA1"))
############ MO plot ###############

MOggplotA1=ggplot(Motif.occurence.allA1,aes(tot.bs,100*Motif.occurence.mean,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA1,Which=="MACPET"),aes(tot.bs,100*Motif.occurence.min),color="red",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA1,Which=="MACPET"),aes(tot.bs,100*Motif.occurence.max),color="red",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA1,Which=="MACS"),aes(tot.bs,100*Motif.occurence.min),color="blue",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA1,Which=="MACS"),aes(tot.bs,100*Motif.occurence.max),color="blue",linetype="dashed",size=0.3)+
    ylab("Motif Occurance (%)")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.x=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0, 0))+ggplot2::scale_y_continuous(expand = c(0, 0),breaks=seq(from=0,to=100,by=5),limits=c(55,85))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.70),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors


MOggplotA1
save("MOggplotA1",file=file.path(DIRsaveA1,"MOggplotA1"))

############ SE plot ###############
SEggplotA1=ggplot(Motif.occurence.allA1,aes(tot.bs,Spatial.error.mean,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA1,Which=="MACPET"),aes(tot.bs,Spatial.error.min),color="red",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA1,Which=="MACPET"),aes(tot.bs,Spatial.error.max),color="red",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA1,Which=="MACS"),aes(tot.bs,Spatial.error.min),color="blue",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA1,Which=="MACS"),aes(tot.bs,Spatial.error.max),color="blue",linetype="dashed",size=0.3)+
   ylab("Spatial Resolution (bp)")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+xlab(" ")+
    ggplot2::scale_x_continuous(expand = c(0, 0))+ggplot2::scale_y_continuous(expand = c(0, 0),breaks=seq(from=0,to=100,by=5),limits=c(15,35))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.80),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

SEggplotA1
save("SEggplotA1",file=file.path(DIRsaveA1,"SEggplotA1"))

############ FDR plot ###############
MACPET_run_1=subset(data.MACPET.Allruns,Run==1)
MACPET_run_1$FDR=MACPET_run_1$FDR*100/max(MACPET_run_1$FDR)
MACS_run_1=subset(data.MACS.Allruns,Run==1)
MACS_run_1$p.values.FDR=MACS_run_1$p.values.FDR*100/max(MACS_run_1$p.values.FDR)
data.FDRA1=data.frame(FDR=c(sort(MACPET_run_1$FDR),sort(MACS_run_1$p.values.FDR)),
                    Which=c(rep("MACPET",nrow(MACPET_run_1)),rep("MACS",nrow(MACS_run_1))),bs=c(1:nrow(MACPET_run_1),1:nrow(MACS_run_1)))
save("data.FDRA1",file=file.path(DIRsaveA1,"data.FDRA1"))

FDRggplotA1=ggplot(data.FDRA1,aes(bs,FDR,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    ylab("FDR (%)")+xlab(" ")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    ggplot2::scale_x_continuous(expand = c(0, 0))+ggplot2::scale_y_continuous(expand = c(0.005, 0),breaks=seq(from=0,to=100,by=10),limits=c(0,100))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.55,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

FDRggplotA1
save("FDRggplotA1",file=file.path(DIRsaveA1,"FDRggplotA1"))

#----------------
#----------------A2: ENCSR000CAD
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
############ Load data MACPET ###############
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAD"
Repetition="REP1"
DIRsaveA2=file.path(SaveDir,"ENCSR000CAD")
#save in the external drive:
DIR=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results",sep="")# the data directory
setwd(DIR)
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACPET.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACPET")
    data.MACPET.run_i=data.motif.overlap.expected
    data.MACPET.run_i$Run=i
    data.MACPET.Allruns=rbind(data.MACPET.Allruns,data.MACPET.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACPET.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACPET.run_1=subset(data.MACPET.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACPET.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACPET.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACPET.run_2=subset(data.MACPET.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACPET.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACPET.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACPET.run_3=subset(data.MACPET.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACPET.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACPET.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACPET.run_4=subset(data.MACPET.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACPET.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACPET.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACPET.run_5=subset(data.MACPET.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACPET.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACPET.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACPET.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                      Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                      Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                      Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                      Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                      Spatial.error.run_1=Spatial.error.run_1,
                                                      Spatial.error.run_2=Spatial.error.run_2,
                                                      Spatial.error.run_3=Spatial.error.run_3,
                                                      Spatial.error.run_4=Spatial.error.run_4,
                                                      Spatial.error.run_5=Spatial.error.run_5,
                                                      Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                      Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                      Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                  tot.motifs.run_2/tot.bs,
                                                                                  tot.motifs.run_3/tot.bs,
                                                                                  tot.motifs.run_4/tot.bs,
                                                                                  tot.motifs.run_5/tot.bs)),
                                                      Spatial.error.min=min(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                      Spatial.error.max=max(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                      Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                                Spatial.error.run_2,
                                                                                Spatial.error.run_3,
                                                                                Spatial.error.run_4,
                                                                                Spatial.error.run_5)),
                                                      tot.bs=x,Which="MACPET")

    # merge:
    data.MACPET.occurrence.spatial.error=rbind(data.MACPET.occurrence.spatial.error,
                                               data.MACPET.occurrence.spatial.error.x)

}
setwd("../../../..")
getwd()
############ Load data MACS ###############
setwd("MACS_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS.run_i=data.motif.overlap.expected
    data.MACS.run_i$Run=i
    data.MACS.Allruns=rbind(data.MACS.Allruns,data.MACS.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                    Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                    Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                    Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                    Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                    Spatial.error.run_1=Spatial.error.run_1,
                                                    Spatial.error.run_2=Spatial.error.run_2,
                                                    Spatial.error.run_3=Spatial.error.run_3,
                                                    Spatial.error.run_4=Spatial.error.run_4,
                                                    Spatial.error.run_5=Spatial.error.run_5,
                                                    Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                              tot.motifs.run_2/tot.bs,
                                                                              tot.motifs.run_3/tot.bs,
                                                                              tot.motifs.run_4/tot.bs,
                                                                              tot.motifs.run_5/tot.bs)),
                                                    Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                              tot.motifs.run_2/tot.bs,
                                                                              tot.motifs.run_3/tot.bs,
                                                                              tot.motifs.run_4/tot.bs,
                                                                              tot.motifs.run_5/tot.bs)),
                                                    Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                    Spatial.error.min=min(c(Spatial.error.run_1,
                                                                            Spatial.error.run_2,
                                                                            Spatial.error.run_3,
                                                                            Spatial.error.run_4,
                                                                            Spatial.error.run_5)),
                                                    Spatial.error.max=max(c(Spatial.error.run_1,
                                                                            Spatial.error.run_2,
                                                                            Spatial.error.run_3,
                                                                            Spatial.error.run_4,
                                                                            Spatial.error.run_5)),
                                                    Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                    tot.bs=x,Which="MACS")

    # merge:
    data.MACS.occurrence.spatial.error=rbind(data.MACS.occurrence.spatial.error,
                                             data.MACS.occurrence.spatial.error.x)

}
setwd("../../..")
getwd()
# setwd("MotifComparizons")
############ Create the whole data ###############
Motif.occurence.allA2=rbind(data.MACPET.occurrence.spatial.error,data.MACS.occurrence.spatial.error)
save("Motif.occurence.allA2",file=file.path(DIRsaveA2,"Motif.occurence.allA2"))

############ MO plot ###############
MOggplotA2=ggplot(Motif.occurence.allA2,aes(tot.bs,100*Motif.occurence.mean,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA2,Which=="MACPET"),aes(tot.bs,100*Motif.occurence.min),color="red",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA2,Which=="MACPET"),aes(tot.bs,100*Motif.occurence.max),color="red",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA2,Which=="MACS"),aes(tot.bs,100*Motif.occurence.min),color="blue",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA2,Which=="MACS"),aes(tot.bs,100*Motif.occurence.max),color="blue",linetype="dashed",size=0.3)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.y=element_blank(),
          axis.title.x=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0, 0),breaks=seq(from=0,to=100,by=3),limits=c(80,95))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

MOggplotA2
save("MOggplotA2",file=file.path(DIRsaveA2,"MOggplotA2"))

############ SE plot ###############
SEggplotA2=ggplot(Motif.occurence.allA2,aes(tot.bs,Spatial.error.mean,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA2,Which=="MACPET"),aes(tot.bs,Spatial.error.min),color="red",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA2,Which=="MACPET"),aes(tot.bs,Spatial.error.max),color="red",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA2,Which=="MACS"),aes(tot.bs,Spatial.error.min),color="blue",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA2,Which=="MACS"),aes(tot.bs,Spatial.error.max),color="blue",linetype="dashed",size=0.3)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.y=element_blank())+xlab("Sorted Binding Sites")+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0, 0),breaks=seq(from=0,to=100,by=3),limits=c(10,30))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.70),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

SEggplotA2
save("SEggplotA2",file=file.path(DIRsaveA2,"SEggplotA2"))

############ FDR plot ###############
MACPET_run_1=subset(data.MACPET.Allruns,Run==1)
MACPET_run_1$FDR=MACPET_run_1$FDR*100/max(MACPET_run_1$FDR)
MACS_run_1=subset(data.MACS.Allruns,Run==1)
MACS_run_1$p.values.FDR=MACS_run_1$p.values.FDR*100/max(MACS_run_1$p.values.FDR)
data.FDRA2=data.frame(FDR=c(sort(MACPET_run_1$FDR),sort(MACS_run_1$p.values.FDR)),
                      Which=c(rep("MACPET",nrow(MACPET_run_1)),rep("MACS",nrow(MACS_run_1))),bs=c(1:nrow(MACPET_run_1),1:nrow(MACS_run_1)))
save("data.FDRA2",file=file.path(DIRsaveA2,"data.FDRA2"))

FDRggplotA2=ggplot(data.FDRA2,aes(bs,FDR,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    xlab("Sorted Binding Sites")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.y=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0.005, 0),breaks=seq(from=0,to=100,by=10),limits=c(0,100))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.55,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

FDRggplotA2
save("FDRggplotA2",file=file.path(DIRsaveA2,"FDRggplotA2"))

#----------------
#----------------A3: ENCSR000CAC
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
############ Load data MACPET ###############
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAC"
Repetition="REP1"
DIRsaveA3=file.path(SaveDir,"ENCSR000CAC")
#save in the external drive:
DIR=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results",sep="")# the data directory
setwd(DIR)
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACPET.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACPET")
    data.MACPET.run_i=data.motif.overlap.expected
    data.MACPET.run_i$Run=i
    data.MACPET.Allruns=rbind(data.MACPET.Allruns,data.MACPET.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACPET.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACPET.run_1=subset(data.MACPET.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACPET.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACPET.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACPET.run_2=subset(data.MACPET.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACPET.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACPET.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACPET.run_3=subset(data.MACPET.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACPET.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACPET.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACPET.run_4=subset(data.MACPET.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACPET.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACPET.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACPET.run_5=subset(data.MACPET.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACPET.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACPET.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACPET.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                      Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                      Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                      Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                      Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                      Spatial.error.run_1=Spatial.error.run_1,
                                                      Spatial.error.run_2=Spatial.error.run_2,
                                                      Spatial.error.run_3=Spatial.error.run_3,
                                                      Spatial.error.run_4=Spatial.error.run_4,
                                                      Spatial.error.run_5=Spatial.error.run_5,
                                                      Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                      Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                      Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                  tot.motifs.run_2/tot.bs,
                                                                                  tot.motifs.run_3/tot.bs,
                                                                                  tot.motifs.run_4/tot.bs,
                                                                                  tot.motifs.run_5/tot.bs)),
                                                      Spatial.error.min=min(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                      Spatial.error.max=max(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                      Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                                Spatial.error.run_2,
                                                                                Spatial.error.run_3,
                                                                                Spatial.error.run_4,
                                                                                Spatial.error.run_5)),
                                                      tot.bs=x,Which="MACPET")

    # merge:
    data.MACPET.occurrence.spatial.error=rbind(data.MACPET.occurrence.spatial.error,
                                               data.MACPET.occurrence.spatial.error.x)

}
setwd("../../../..")
getwd()
############ Load data MACS ###############
setwd("MACS_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS.run_i=data.motif.overlap.expected
    data.MACS.run_i$Run=i
    data.MACS.Allruns=rbind(data.MACS.Allruns,data.MACS.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                    Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                    Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                    Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                    Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                    Spatial.error.run_1=Spatial.error.run_1,
                                                    Spatial.error.run_2=Spatial.error.run_2,
                                                    Spatial.error.run_3=Spatial.error.run_3,
                                                    Spatial.error.run_4=Spatial.error.run_4,
                                                    Spatial.error.run_5=Spatial.error.run_5,
                                                    Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                              tot.motifs.run_2/tot.bs,
                                                                              tot.motifs.run_3/tot.bs,
                                                                              tot.motifs.run_4/tot.bs,
                                                                              tot.motifs.run_5/tot.bs)),
                                                    Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                              tot.motifs.run_2/tot.bs,
                                                                              tot.motifs.run_3/tot.bs,
                                                                              tot.motifs.run_4/tot.bs,
                                                                              tot.motifs.run_5/tot.bs)),
                                                    Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                    Spatial.error.min=min(c(Spatial.error.run_1,
                                                                            Spatial.error.run_2,
                                                                            Spatial.error.run_3,
                                                                            Spatial.error.run_4,
                                                                            Spatial.error.run_5)),
                                                    Spatial.error.max=max(c(Spatial.error.run_1,
                                                                            Spatial.error.run_2,
                                                                            Spatial.error.run_3,
                                                                            Spatial.error.run_4,
                                                                            Spatial.error.run_5)),
                                                    Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                    tot.bs=x,Which="MACS")

    # merge:
    data.MACS.occurrence.spatial.error=rbind(data.MACS.occurrence.spatial.error,
                                             data.MACS.occurrence.spatial.error.x)

}
setwd("../../..")
getwd()
# setwd("MotifComparizons")
############ Create the whole data ###############
Motif.occurence.allA3=rbind(data.MACPET.occurrence.spatial.error,data.MACS.occurrence.spatial.error)
save("Motif.occurence.allA3",file=file.path(DIRsaveA3,"Motif.occurence.allA3"))

############ MO plot ###############
MOggplotA3=ggplot(Motif.occurence.allA3,aes(tot.bs,100*Motif.occurence.mean,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA3,Which=="MACPET"),aes(tot.bs,100*Motif.occurence.min),color="red",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA3,Which=="MACPET"),aes(tot.bs,100*Motif.occurence.max),color="red",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA3,Which=="MACS"),aes(tot.bs,100*Motif.occurence.min),color="blue",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA3,Which=="MACS"),aes(tot.bs,100*Motif.occurence.max),color="blue",linetype="dashed",size=0.3)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.y=element_blank(),
          axis.title.x=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0, 0),breaks=seq(from=0,to=100,by=1),limits=c(92,98))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.70),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

MOggplotA3
save("MOggplotA3",file=file.path(DIRsaveA3,"MOggplotA3"))

############ SE plot ###############
SEggplotA3=ggplot(Motif.occurence.allA3,aes(tot.bs,Spatial.error.mean,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA3,Which=="MACPET"),aes(tot.bs,Spatial.error.min),color="red",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA3,Which=="MACPET"),aes(tot.bs,Spatial.error.max),color="red",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA3,Which=="MACS"),aes(tot.bs,Spatial.error.min),color="blue",linetype="dashed",size=0.3)+
    # geom_line(data=subset(Motif.occurence.allA3,Which=="MACS"),aes(tot.bs,Spatial.error.max),color="blue",linetype="dashed",size=0.3)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.y=element_blank())+xlab(" ")+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0, 0),breaks=seq(from=0,to=100,by=2),limits=c(10,18))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.70),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

SEggplotA3
save("SEggplotA3",file=file.path(DIRsaveA3,"SEggplotA3"))

############ FDR plot ###############
MACPET_run_1=subset(data.MACPET.Allruns,Run==1)
MACPET_run_1$FDR=MACPET_run_1$FDR*100/max(MACPET_run_1$FDR)
MACS_run_1=subset(data.MACS.Allruns,Run==1)
MACS_run_1$p.values.FDR=MACS_run_1$p.values.FDR*100/max(MACS_run_1$p.values.FDR)
data.FDRA3=data.frame(FDR=c(sort(MACPET_run_1$FDR),sort(MACS_run_1$p.values.FDR)),
                      Which=c(rep("MACPET",nrow(MACPET_run_1)),rep("MACS",nrow(MACS_run_1))),bs=c(1:nrow(MACPET_run_1),1:nrow(MACS_run_1)))
save("data.FDRA3",file=file.path(DIRsaveA3,"data.FDRA3"))

FDRggplotA3=ggplot(data.FDRA3,aes(bs,FDR,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    xlab("")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.y=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0.005, 0),breaks=seq(from=0,to=100,by=10),limits=c(0,100))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.55,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

FDRggplotA3
save("FDRggplotA3",file=file.path(DIRsaveA3,"FDRggplotA3"))

#################################################################################################################
################################ Figures for FDR for the rest: ######################################
#################################################################################################################
# The saving directory
#----------------
#----------------A4: ENCSR000FDD
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000FDD"
Repetition="REP1"
DIRsaveA4=file.path(SaveDir,"ENCSR000FDD")
############ Load data MACPET ###############
DIRMACPET=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results/MACPET_psfitData",sep="")# the data directory
load(DIRMACPET)
PeaksMACPET=metadata(MACPET_psfitData)
PeaksMACPET=PeaksMACPET$Peaks.Info
PeaksMACPET=PeaksMACPET[with(PeaksMACPET,order(FDR,decreasing=F)),]
PeaksMACPET=PeaksMACPET[1:5000,]
PeaksMACPET$FDR=PeaksMACPET$FDR*100/max(PeaksMACPET$FDR)
############ Load data MACPET ###############
DIRMACS=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACS_results/MACS.peaks.out_peaks.update",sep="")# the data directory
load(DIRMACS)
PeaksMACS=MACS.peaks.out_peaks.update[with(MACS.peaks.out_peaks.update,order(p.values.FDR,decreasing=F)),]
PeaksMACS=PeaksMACS[1:5000,]
PeaksMACS$p.values.FDR=PeaksMACS$p.values.FDR*100/max(PeaksMACS$p.values.FDR)
# plot:
data.FDRA4=data.frame(FDR=c(sort(PeaksMACPET$FDR),sort(PeaksMACS$p.values.FDR)),
                      Which=c(rep("MACPET",nrow(PeaksMACPET)),rep("MACS",nrow(PeaksMACS))),bs=c(1:nrow(PeaksMACPET),1:nrow(PeaksMACS)))
save("data.FDRA4",file=file.path(DIRsaveA4,"data.FDRA4"))

FDRggplotA4=ggplot(data.FDRA4,aes(bs,FDR,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    xlab("")+ylab("FDR (%)")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0.005, 0),breaks=seq(from=0,to=100,by=10),limits=c(0,100))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.55,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

FDRggplotA4
save("FDRggplotA4",file=file.path(DIRsaveA4,"FDRggplotA4"))

#----------------
#----------------A5: ENCSR000FDG or Not_specified_1
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000FDG"
Repetition="REP1"
DIRsaveA5=file.path(SaveDir,"ENCSR000FDG")
############ Load data MACPET ###############
DIRMACPET=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results/MACPET_psfitData",sep="")# the data directory
load(DIRMACPET)
PeaksMACPET=metadata(MACPET_psfitData)
PeaksMACPET=PeaksMACPET$Peaks.Info
PeaksMACPET=PeaksMACPET[with(PeaksMACPET,order(FDR,decreasing=F)),]
PeaksMACPET=PeaksMACPET[1:5000,]
PeaksMACPET$FDR=PeaksMACPET$FDR*100/max(PeaksMACPET$FDR)
############ Load data MACPET ###############
DIRMACS=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACS_results/MACS.peaks.out_peaks.update",sep="")# the data directory
load(DIRMACS)
PeaksMACS=MACS.peaks.out_peaks.update[with(MACS.peaks.out_peaks.update,order(p.values.FDR,decreasing=F)),]
PeaksMACS=PeaksMACS[1:5000,]
PeaksMACS$p.values.FDR=PeaksMACS$p.values.FDR*100/max(PeaksMACS$p.values.FDR)
# plot:
data.FDRA5=data.frame(FDR=c(sort(PeaksMACPET$FDR),sort(PeaksMACS$p.values.FDR)),
                      Which=c(rep("MACPET",nrow(PeaksMACPET)),rep("MACS",nrow(PeaksMACS))),bs=c(1:nrow(PeaksMACPET),1:nrow(PeaksMACS)))
save("data.FDRA5",file=file.path(DIRsaveA5,"data.FDRA5"))

FDRggplotA5=ggplot(data.FDRA5,aes(bs,FDR,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    xlab("Sorted Binding Sites")+ylab("")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.y=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0.005, 0),breaks=seq(from=0,to=100,by=10),limits=c(0,100))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.55,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

FDRggplotA5
save("FDRggplotA5",file=file.path(DIRsaveA5,"FDRggplotA5"))

#----------------
#----------------A6: ENCSR000BZY
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZY"
Repetition="REP1"
DIRsaveA6=file.path(SaveDir,"ENCSR000BZY")
############ Load data MACPET ###############
DIRMACPET=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results/MACPET_psfitData",sep="")# the data directory
load(DIRMACPET)
PeaksMACPET=metadata(MACPET_psfitData)
PeaksMACPET=PeaksMACPET$Peaks.Info
PeaksMACPET=PeaksMACPET[with(PeaksMACPET,order(FDR,decreasing=F)),]
PeaksMACPET=PeaksMACPET[1:5000,]
PeaksMACPET$FDR=PeaksMACPET$FDR*100/max(PeaksMACPET$FDR)
############ Load data MACPET ###############
DIRMACS=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACS_results/MACS.peaks.out_peaks.update",sep="")# the data directory
load(DIRMACS)
PeaksMACS=MACS.peaks.out_peaks.update[with(MACS.peaks.out_peaks.update,order(p.values.FDR,decreasing=F)),]
PeaksMACS=PeaksMACS[1:5000,]
PeaksMACS$p.values.FDR=PeaksMACS$p.values.FDR*100/max(PeaksMACS$p.values.FDR)
# plot:
data.FDRA6=data.frame(FDR=c(sort(PeaksMACPET$FDR),sort(PeaksMACS$p.values.FDR)),
                      Which=c(rep("MACPET",nrow(PeaksMACPET)),rep("MACS",nrow(PeaksMACS))),bs=c(1:nrow(PeaksMACPET),1:nrow(PeaksMACS)))
save("data.FDRA6",file=file.path(DIRsaveA6,"data.FDRA6"))

FDRggplotA6=ggplot(data.FDRA6,aes(bs,FDR,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    xlab("")+ylab("")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.y=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0.005, 0),breaks=seq(from=0,to=100,by=10),limits=c(0,100))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.55,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

FDRggplotA6
save("FDRggplotA6",file=file.path(DIRsaveA6,"FDRggplotA6"))


#################################################################################################################
################################ Figures for Interactions MANGO: ###############################
#################################################################################################################
#-------------------------------------------------------
#----------------A1: ENCSR000BZZ
#-------------------------------------------------------
# The saving directory
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZZ"
Repetition="REP1"
DIRsaveA1=file.path(SaveDir,"ENCSR000BZZ/MANGO_interactions")
if(!dir.exists(DIRsaveA1)) dir.create(DIRsaveA1)
DIRdata=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/InteractionsResults_MANGO_MACPETInitials/",sep="")# the data directory
############ Load data MACPET ###############
MACPET_files=list.files(file.path(DIRdata,"MACPET"),pattern="Interactions_")
MACPET_Sign_Int_MANGO_data_A1=data.frame()#for significant interactions
MACPET_HeatMap_mat_MANGO_A1=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACPET_HeatMap_mat_MANGO_A1)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACPET_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACPET_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACPET",folder,"MACPET_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACPET_slopPeak_folder=MACPET_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACPET_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACPET_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACPET_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACPET_slopPeak_folder$chr,
                                                  ranges=IRanges::IRanges(start=MACPET_slopPeak_folder$start,
                                                                          end=MACPET_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACPET",folder,"MACPET.interactions.fdr.mango"))){
        MACPET_Interactions_folder=utils::read.table(file.path(DIRdata,"MACPET",folder,"MACPET.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACPET_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACPET_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACPET_Interactions_folder_update=data.frame(chr=c(MACPET_Interactions_folder$chr_from,
                                                           MACPET_Interactions_folder$chr_to),
                                                     start=c(MACPET_Interactions_folder$start_from,
                                                             MACPET_Interactions_folder$start_to),
                                                     end=c(MACPET_Interactions_folder$end_from,
                                                           MACPET_Interactions_folder$end_to))
        # make granges
        MACPET_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACPET_Interactions_folder_update$chr,
                                                                 ranges=IRanges::IRanges(start=MACPET_Interactions_folder_update$start,
                                                                                         end=MACPET_Interactions_folder_update$end))

        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACPET_slopPeak_folder,MACPET_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACPET_HeatMap_mat_MANGO_A1)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACPET_HeatMap_mat_MANGO_A1[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACPET",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#else no interacions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACPET",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }

    MACPET_Sign_Int_MANGO_data_A1=rbind(MACPET_Sign_Int_MANGO_data_A1,Interactions_folder)
}
# sort by merging window
MACPET_Sign_Int_MANGO_data_A1=MACPET_Sign_Int_MANGO_data_A1[with(MACPET_Sign_Int_MANGO_data_A1,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
MACPET_HeatMap_mat_MANGO_A1=MACPET_HeatMap_mat_MANGO_A1[c(1:(Top_max_used)),]#remove rows
MACPET_HeatMap_mat_MANGO_A1[which(is.na(MACPET_HeatMap_mat_MANGO_A1),arr.ind=TRUE)]=0
save("MACPET_HeatMap_mat_MANGO_A1",file=file.path(DIRsaveA1,"MACPET_HeatMap_mat_MANGO_A1"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
MACPETInteractions500_MANGO_A1=utils::read.table(file.path(DIRdata,"MACPET","Interactions_500","MACPET.interactions.fdr.mango"),stringsAsFactors=FALSE)
colnames(MACPETInteractions500_MANGO_A1)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
save("MACPETInteractions500_MANGO_A1",file=file.path(DIRsaveA1,"MACPETInteractions500_MANGO_A1"))

############ Load data MACS ###############
MACS_files=list.files(file.path(DIRdata,"MACS"),pattern="Interactions_")
MACS_0.05_Sign_Int_MANGO_data_A1=data.frame()
MACS_HeatMap_mat_MANGO_A1=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A1)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A1)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A1[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }

        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.05",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.05",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }

    MACS_0.05_Sign_Int_MANGO_data_A1=rbind(MACS_0.05_Sign_Int_MANGO_data_A1,Interactions_folder)
}
# sort by merging window
MACS_0.05_Sign_Int_MANGO_data_A1=MACS_0.05_Sign_Int_MANGO_data_A1[with(MACS_0.05_Sign_Int_MANGO_data_A1,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
MACS_HeatMap_mat_MANGO_A1=MACS_HeatMap_mat_MANGO_A1[c(1:(Top_max_used)),]#remove rows
MACS_HeatMap_mat_MANGO_A1[which(is.na(MACS_HeatMap_mat_MANGO_A1),arr.ind=TRUE)]=0
save("MACS_HeatMap_mat_MANGO_A1",file=file.path(DIRsaveA1,"MACS_HeatMap_mat_MANGO_A1"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
MACSInteractions500_MANGO_A1=utils::read.table(file.path(DIRdata,"MACS","Interactions_500","MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
colnames(MACSInteractions500_MANGO_A1)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
save("MACSInteractions500_MANGO_A1",file=file.path(DIRsaveA1,"MACSInteractions500_MANGO_A1"))

############ Load data MACS 0.01 ###############
MACS_files_0.01=list.files(file.path(DIRdata,"MACS_0.01"),pattern="Interactions_")
MACS_0.01_Sign_Int_MANGO_data_A1=data.frame()
MACS_HeatMap_mat_MANGO_A1=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A1)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_files_0.01){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS_0.01",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS_0.01",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS_0.01",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A1)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A1[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.01",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.01",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }
    
    MACS_0.01_Sign_Int_MANGO_data_A1=rbind(MACS_0.01_Sign_Int_MANGO_data_A1,Interactions_folder)
}
# sort by merging window
MACS_0.01_Sign_Int_MANGO_data_A1=MACS_0.01_Sign_Int_MANGO_data_A1[with(MACS_0.01_Sign_Int_MANGO_data_A1,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
MACS_HeatMap_mat_MANGO_A1=MACS_HeatMap_mat_MANGO_A1[c(1:(Top_max_used)),]#remove rows
MACS_HeatMap_mat_MANGO_A1[which(is.na(MACS_HeatMap_mat_MANGO_A1),arr.ind=TRUE)]=0
save("MACS_HeatMap_mat_MANGO_A1",file=file.path(DIRsaveA1,"MACS_HeatMap_mat_MANGO_A1"))
# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
MACSInteractions500_MANGO_A1=utils::read.table(file.path(DIRdata,"MACS","Interactions_500","MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
colnames(MACSInteractions500_MANGO_A1)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
save("MACSInteractions500_MANGO_A1",file=file.path(DIRsaveA1,"MACSInteractions500_MANGO_A1"))

############ Load data MACS N ###############
MACS_files_N=list.files(file.path(DIRdata,"MACS_N"),pattern="Interactions_")
MACS_N_Sign_Int_MANGO_data_A1=data.frame()
MACS_HeatMap_mat_MANGO_A1=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A1)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_files_N){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS_N",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS_N",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS_N",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A1)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A1[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS N",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS N",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }
    
    MACS_N_Sign_Int_MANGO_data_A1=rbind(MACS_N_Sign_Int_MANGO_data_A1,Interactions_folder)
}
# sort by merging window
MACS_N_Sign_Int_MANGO_data_A1=MACS_N_Sign_Int_MANGO_data_A1[with(MACS_N_Sign_Int_MANGO_data_A1,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
MACS_HeatMap_mat_MANGO_A1=MACS_HeatMap_mat_MANGO_A1[c(1:(Top_max_used)),]#remove rows
MACS_HeatMap_mat_MANGO_A1[which(is.na(MACS_HeatMap_mat_MANGO_A1),arr.ind=TRUE)]=0
save("MACS_HeatMap_mat_MANGO_A1",file=file.path(DIRsaveA1,"MACS_HeatMap_mat_MANGO_A1"))
# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
MACSInteractions500_MANGO_A1=utils::read.table(file.path(DIRdata,"MACS","Interactions_500","MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
colnames(MACSInteractions500_MANGO_A1)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
save("MACSInteractions500_MANGO_A1",file=file.path(DIRsaveA1,"MACSInteractions500_MANGO_A1"))



############ Merge and make plot ###############
MACPET_MACS_Sign_Int_MANGO_data_A1=rbind(MACPET_Sign_Int_MANGO_data_A1,MACS_0.05_Sign_Int_MANGO_data_A1,MACS_0.01_Sign_Int_MANGO_data_A1,MACS_N_Sign_Int_MANGO_data_A1)
save("MACPET_MACS_Sign_Int_MANGO_data_A1",file=file.path(DIRsaveA1,"MACPET_MACS_Sign_Int_MANGO_data_A1"))

# ------------
# significant interactions plot (in row 1):
# ------------
SigInteractions_MANGO_plot_A1=ggplot(MACPET_MACS_Sign_Int_MANGO_data_A1,aes(Merge_Window,Total_significant_int,color=factor(Which),linetype=factor(Which)))+
    geom_line(size=0.3)+ylab("Significant Interactions")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.x=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.005, 0),breaks=seq(0,1000,100))+
    ggplot2::scale_y_continuous(expand = c(0.01,0),breaks=seq(min(MACPET_MACS_Sign_Int_MANGO_data_A1$Total_significant_int),max(MACPET_MACS_Sign_Int_MANGO_data_A1$Total_significant_int),25))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(1.15,0.3),legend.title=element_blank(),
          legend.text=element_text(size=5))+guides(col = guide_legend(nrow=2))+
    scale_colour_manual(values=c("red","blue","green","brown"))#coloring factors

SigInteractions_MANGO_plot_A1
save("SigInteractions_MANGO_plot_A1",file=file.path(DIRsaveA1,"SigInteractions_MANGO_plot_A1"))


# ------------
# plots for percentage of peaks inn and peaks used in the interactions(in row 2)
# ------------
PeaksPercentegeInvolved_MANGO_plot_A1=ggplot(MACPET_MACS_Sign_Int_MANGO_data_A1,aes(Merge_Window,100*Peaks_involved_1/Peaks_in,color=factor(Which),linetype=factor(Which)))+
    geom_line(size=0.3)+ylab("Peaks used (%)")+xlab(" ")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    ggplot2::scale_x_continuous(expand = c(0.005, 0),breaks=seq(0,1000,100))+
    ggplot2::scale_y_continuous(expand = c(0.01,0),breaks=seq(min(round(100*MACPET_MACS_Sign_Int_MANGO_data_A1$Peaks_involved_1/MACPET_MACS_Sign_Int_MANGO_data_A1$Peaks_in)),
                                                           max(round(100*MACPET_MACS_Sign_Int_MANGO_data_A1$Peaks_involved_1/MACPET_MACS_Sign_Int_MANGO_data_A1$Peaks_in)),0.50))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.35),legend.title=element_blank(),
          legend.text=element_text(size=6))+guides(col = guide_legend(nrow=2))+
    scale_colour_manual(values=c("red","blue","green","brown"))#coloring factors

PeaksPercentegeInvolved_MANGO_plot_A1
save("PeaksPercentegeInvolved_MANGO_plot_A1",file=file.path(DIRsaveA1,"PeaksPercentegeInvolved_MANGO_plot_A1"))

# ------------
# Plot for interactions overlap for 500 window:
# ------------
# find interactions overlap:
# make Granges MACPET:
Anchor1_MACPET_MANGO_A1=GenomicRanges::GRanges(seqnames=MACPETInteractions500_MANGO_A1$chr_from,
                                         ranges=IRanges::IRanges(start=MACPETInteractions500_MANGO_A1$start_from,
                                                                 end=MACPETInteractions500_MANGO_A1$end_from))
Anchor2_MACPET_MANGO_A1=GenomicRanges::GRanges(seqnames=MACPETInteractions500_MANGO_A1$chr_to,
                                         ranges=IRanges::IRanges(start=MACPETInteractions500_MANGO_A1$start_to,
                                                                 end=MACPETInteractions500_MANGO_A1$end_to))
# make interactions
GInteractions_MACPET_MANGO_A1=GInteractions(Anchor1_MACPET_MANGO_A1,Anchor2_MACPET_MANGO_A1)
# make Granges MACS:
Anchor1_MACS_MANGO_A1=GenomicRanges::GRanges(seqnames=MACSInteractions500_MANGO_A1$chr_from,
                                       ranges=IRanges::IRanges(start=MACSInteractions500_MANGO_A1$start_from,
                                                               end=MACSInteractions500_MANGO_A1$end_from))
Anchor2_MACS_MANGO_A1=GenomicRanges::GRanges(seqnames=MACSInteractions500_MANGO_A1$chr_to,
                                       ranges=IRanges::IRanges(start=MACSInteractions500_MANGO_A1$start_to,
                                                               end=MACSInteractions500_MANGO_A1$end_to))

GInteractions_MACS_MANGO_A1=GInteractions(Anchor1_MACS_MANGO_A1,Anchor2_MACS_MANGO_A1)


# make calculations for creating the venn data:
library(ggforce)
TotInteractions_MACPET_MANGO_A1=length(GInteractions_MACPET_MANGO_A1)
TotInteractions_MACS_MANGO_A1=length(GInteractions_MACS_MANGO_A1)
TotInteractions_MANGO_A1=TotInteractions_MACPET_MANGO_A1+TotInteractions_MACS_MANGO_A1
TotCommonInteractions_MACStoMACPET_MANGO_A1=length(InteractionSet::findOverlaps(GInteractions_MACPET_MANGO_A1,GInteractions_MACS_MANGO_A1))
# change to 100 scale
TotInteractions_MACPET_MANGO_A1_100=TotInteractions_MACPET_MANGO_A1*100/TotInteractions_MANGO_A1
TotInteractions_MACS_MANGO_A1_100=TotInteractions_MACS_MANGO_A1*100/TotInteractions_MANGO_A1
TotCommonInteractions_MACStoMACPET_MANGO_A1_100=TotCommonInteractions_MACStoMACPET_MANGO_A1*100/TotInteractions_MANGO_A1

VennInteractions_MANGO_A1_data=data.frame(x0=c(0+TotInteractions_MACPET_MANGO_A1_100/2+TotCommonInteractions_MACStoMACPET_MANGO_A1_100/2,#macpet
                                 100-TotInteractions_MACS_MANGO_A1_100/2-TotCommonInteractions_MACStoMACPET_MANGO_A1_100/2,NA),#macs
                            y0=c(0.5,0.5,NA),
                            r=c(TotInteractions_MACPET_MANGO_A1_100/2,TotInteractions_MACS_MANGO_A1_100/2,NA),
                            Totals=c(TotInteractions_MACPET_MANGO_A1,TotInteractions_MACS_MANGO_A1,TotCommonInteractions_MACStoMACPET_MANGO_A1),
                            color=c("red","blue",NA),which=c("MACPET","MACS",NA))

save("VennInteractions_MANGO_A1_data",file=file.path(DIRsaveA1,"VennInteractions_MANGO_A1_data"))


VennInteractions_Overlap_MANGO_A1_plot=ggplot()+geom_circle(data=VennInteractions_MANGO_A1_data[2,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="blue",alpha=0.6,size=0.2)+
    geom_circle(data=VennInteractions_MANGO_A1_data[1,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="red",alpha=0.6,size=0.2)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),axis.ticks=element_blank(),
          axis.text=element_blank(),axis.title=element_blank())+
    # add numbers for overlaps:
    annotate("text", x=50, y=0.5, label= as.character(VennInteractions_MANGO_A1_data$Totals[3]),size=2.5)+#common
    annotate("text", x=100-VennInteractions_MANGO_A1_data$r[2], y=0.5, label= as.character(VennInteractions_MANGO_A1_data$Totals[2]-VennInteractions_MANGO_A1_data$Totals[3]),size=2.5)+#macs
    annotate("text", x=0+VennInteractions_MANGO_A1_data$r[1], y=0.5, label= as.character(VennInteractions_MANGO_A1_data$Totals[1]-VennInteractions_MANGO_A1_data$Totals[3]),size=2.5)+#MACPET
    # add groups:
    annotate("text", x=30, y=30, label= "MACPET",size=3)+#MACPET
    annotate("text", x=70, y=30, label= "MACS",size=3)#MACPET


VennInteractions_Overlap_MANGO_A1_plot
save("VennInteractions_Overlap_MANGO_A1_plot",file=file.path(DIRsaveA1,"VennInteractions_Overlap_MANGO_A1_plot"))

# ------------
# Plot for interactions length for 500 window:
# ------------
scientific_10 <- function(x) {
    y=parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
    if(x[1]==0) y[1]=0

    return(y)
}#for scientific values
# Find interaction Distances:
Interactions_Dist_MACPET_MANGO_A1=pairdist(GInteractions_MACPET_MANGO_A1)
Interactions_Dist_MACS_MANGO_A1=pairdist(GInteractions_MACS_MANGO_A1)
Interaction_Dist_Merged_MANGO_A1=data.frame(Dist=c(Interactions_Dist_MACPET_MANGO_A1,Interactions_Dist_MACS_MANGO_A1),
                                      Which=c(rep("MACPET",length(Interactions_Dist_MACPET_MANGO_A1)),
                                              rep("MACS",length(Interactions_Dist_MACS_MANGO_A1))))
save("Interaction_Dist_Merged_MANGO_A1",file=file.path(DIRsaveA1,"Interaction_Dist_Merged_MANGO_A1"))


Interaction_Dist_MANGO_Plot_A1=ggplot(Interaction_Dist_Merged_MANGO_A1,aes(Dist,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(label=scientific_10,expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    xlab("")+
    ylab("Frequency")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

Interaction_Dist_MANGO_Plot_A1
save("Interaction_Dist_MANGO_Plot_A1",file=file.path(DIRsaveA1,"Interaction_Dist_MANGO_Plot_A1"))


#-------------------------------------------------------
#----------------A2: ENCSR000CAD
#-------------------------------------------------------
# The saving directory
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAD"
Repetition="REP1"
DIRsaveA2=file.path(SaveDir,"ENCSR000CAD/MANGO_interactions")
if(!dir.exists(DIRsaveA2)) dir.create(DIRsaveA2)
DIRdata=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/InteractionsResults_MANGO_MACPETInitials/",sep="")# the data directory
############ Load data MACPET ###############
MACPET_files=list.files(file.path(DIRdata,"MACPET"),pattern="Interactions_")
MACPET_Sign_Int_MANGO_data_A2=data.frame()#for significant interactions
MACPET_HeatMap_mat_MANGO_A2=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACPET_HeatMap_mat_MANGO_A2)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACPET_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACPET_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACPET",folder,"MACPET_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACPET_slopPeak_folder=MACPET_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACPET_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACPET_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACPET_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACPET_slopPeak_folder$chr,
                                                  ranges=IRanges::IRanges(start=MACPET_slopPeak_folder$start,
                                                                          end=MACPET_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACPET",folder,"MACPET.interactions.fdr.mango"))){
        MACPET_Interactions_folder=utils::read.table(file.path(DIRdata,"MACPET",folder,"MACPET.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACPET_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACPET_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACPET_Interactions_folder_update=data.frame(chr=c(MACPET_Interactions_folder$chr_from,
                                                           MACPET_Interactions_folder$chr_to),
                                                     start=c(MACPET_Interactions_folder$start_from,
                                                             MACPET_Interactions_folder$start_to),
                                                     end=c(MACPET_Interactions_folder$end_from,
                                                           MACPET_Interactions_folder$end_to))
        # make granges
        MACPET_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACPET_Interactions_folder_update$chr,
                                                                 ranges=IRanges::IRanges(start=MACPET_Interactions_folder_update$start,
                                                                                         end=MACPET_Interactions_folder_update$end))

        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACPET_slopPeak_folder,MACPET_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACPET_HeatMap_mat_MANGO_A2)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACPET_HeatMap_mat_MANGO_A2[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACPET",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#else no interacions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACPET",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }

    MACPET_Sign_Int_MANGO_data_A2=rbind(MACPET_Sign_Int_MANGO_data_A2,Interactions_folder)
}
# sort by merging window
MACPET_Sign_Int_MANGO_data_A2=MACPET_Sign_Int_MANGO_data_A2[with(MACPET_Sign_Int_MANGO_data_A2,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
MACPET_HeatMap_mat_MANGO_A2=MACPET_HeatMap_mat_MANGO_A2[c(1:(Top_max_used)),]#remove rows
MACPET_HeatMap_mat_MANGO_A2[which(is.na(MACPET_HeatMap_mat_MANGO_A2),arr.ind=TRUE)]=0
save("MACPET_HeatMap_mat_MANGO_A2",file=file.path(DIRsaveA2,"MACPET_HeatMap_mat_MANGO_A2"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
MACPETInteractions500_MANGO_A2=utils::read.table(file.path(DIRdata,"MACPET","Interactions_500","MACPET.interactions.fdr.mango"),stringsAsFactors=FALSE)
colnames(MACPETInteractions500_MANGO_A2)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
save("MACPETInteractions500_MANGO_A2",file=file.path(DIRsaveA2,"MACPETInteractions500_MANGO_A2"))

############ Load data MACS ###############
MACS_files=list.files(file.path(DIRdata,"MACS"),pattern="Interactions_")
MACS_0.05_Sign_Int_MANGO_data_A2=data.frame()
MACS_HeatMap_mat_MANGO_A2=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A2)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A2)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A2[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }

        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.05",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.05",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }

    MACS_0.05_Sign_Int_MANGO_data_A2=rbind(MACS_0.05_Sign_Int_MANGO_data_A2,Interactions_folder)
}
# sort by merging window
MACS_0.05_Sign_Int_MANGO_data_A2=MACS_0.05_Sign_Int_MANGO_data_A2[with(MACS_0.05_Sign_Int_MANGO_data_A2,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
MACS_HeatMap_mat_MANGO_A2=MACS_HeatMap_mat_MANGO_A2[c(1:(Top_max_used)),]#remove rows
MACS_HeatMap_mat_MANGO_A2[which(is.na(MACS_HeatMap_mat_MANGO_A2),arr.ind=TRUE)]=0
save("MACS_HeatMap_mat_MANGO_A2",file=file.path(DIRsaveA2,"MACS_HeatMap_mat_MANGO_A2"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
MACSInteractions500_MANGO_A2=utils::read.table(file.path(DIRdata,"MACS","Interactions_500","MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
colnames(MACSInteractions500_MANGO_A2)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
save("MACSInteractions500_MANGO_A2",file=file.path(DIRsaveA2,"MACSInteractions500_MANGO_A2"))

############ Load data MACS 0.01 ###############
MACS_0.01_files=list.files(file.path(DIRdata,"MACS_0.01"),pattern="Interactions_")
MACS_0.01_Sign_Int_MANGO_data_A2=data.frame()
MACS_HeatMap_mat_MANGO_A2=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A2)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_0.01_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS_0.01",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS_0.01",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS_0.01",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A2)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A2[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.01",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.01",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }
    
    MACS_0.01_Sign_Int_MANGO_data_A2=rbind(MACS_0.01_Sign_Int_MANGO_data_A2,Interactions_folder)
}
# sort by merging window
MACS_0.01_Sign_Int_MANGO_data_A2=MACS_0.01_Sign_Int_MANGO_data_A2[with(MACS_0.01_Sign_Int_MANGO_data_A2,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
# MACS_HeatMap_mat_MANGO_A2=MACS_HeatMap_mat_MANGO_A2[c(1:(Top_max_used)),]#remove rows
# MACS_HeatMap_mat_MANGO_A2[which(is.na(MACS_HeatMap_mat_MANGO_A2),arr.ind=TRUE)]=0
# save("MACS_HeatMap_mat_MANGO_A2",file=file.path(DIRsaveA2,"MACS_HeatMap_mat_MANGO_A2"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
# MACSInteractions500_MANGO_A2=utils::read.table(file.path(DIRdata,"MACS","Interactions_500","MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
# colnames(MACSInteractions500_MANGO_A2)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
# save("MACSInteractions500_MANGO_A2",file=file.path(DIRsaveA2,"MACSInteractions500_MANGO_A2"))

############ Load data MACS N ###############
MACS_N_files=list.files(file.path(DIRdata,"MACS_N"),pattern="Interactions_")
MACS_N_Sign_Int_MANGO_data_A2=data.frame()
MACS_HeatMap_mat_MANGO_A2=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A2)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_N_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS_N",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS_N",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS_N",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A2)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A2[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS N",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS N",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }
    
    MACS_N_Sign_Int_MANGO_data_A2=rbind(MACS_N_Sign_Int_MANGO_data_A2,Interactions_folder)
}
# sort by merging window
MACS_N_Sign_Int_MANGO_data_A2=MACS_N_Sign_Int_MANGO_data_A2[with(MACS_N_Sign_Int_MANGO_data_A2,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
# MACS_HeatMap_mat_MANGO_A2=MACS_HeatMap_mat_MANGO_A2[c(1:(Top_max_used)),]#remove rows
# MACS_HeatMap_mat_MANGO_A2[which(is.na(MACS_HeatMap_mat_MANGO_A2),arr.ind=TRUE)]=0
# save("MACS_HeatMap_mat_MANGO_A2",file=file.path(DIRsaveA2,"MACS_HeatMap_mat_MANGO_A2"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
# MACSInteractions500_MANGO_A2=utils::read.table(file.path(DIRdata,"MACS","Interactions_500","MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
# colnames(MACSInteractions500_MANGO_A2)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
# save("MACSInteractions500_MANGO_A2",file=file.path(DIRsaveA2,"MACSInteractions500_MANGO_A2"))


############ Merge and make plot ###############
MACPET_MACS_Sign_Int_MANGO_data_A2=rbind(MACPET_Sign_Int_MANGO_data_A2,MACS_0.05_Sign_Int_MANGO_data_A2,MACS_0.01_Sign_Int_MANGO_data_A2,MACS_N_Sign_Int_MANGO_data_A2)
save("MACPET_MACS_Sign_Int_MANGO_data_A2",file=file.path(DIRsaveA2,"MACPET_MACS_Sign_Int_MANGO_data_A2"))

# ------------
# significant interactions plot (in row 1):
# ------------
SigInteractions_MANGO_plot_A2=ggplot(MACPET_MACS_Sign_Int_MANGO_data_A2,aes(Merge_Window,Total_significant_int,color=factor(Which),linetype=factor(Which)))+
    geom_line(size=0.3)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.x=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.005, 0),breaks=seq(0,1000,100))+
    ggplot2::scale_y_continuous(expand = c(0.01,0),breaks=seq(min(MACPET_MACS_Sign_Int_MANGO_data_A2$Total_significant_int),max(MACPET_MACS_Sign_Int_MANGO_data_A2$Total_significant_int),800))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(1.05,.05),legend.title=element_blank(),
          legend.text=element_text(size=6),axis.title.x=element_blank(),axis.title.y=element_blank())+guides(col = guide_legend(nrow=2))+
    scale_colour_manual(values=c("red","blue","green","brown"))#coloring factors

SigInteractions_MANGO_plot_A2
save("SigInteractions_MANGO_plot_A2",file=file.path(DIRsaveA2,"SigInteractions_MANGO_plot_A2"))


# ------------
# plots for percentage of peaks inn and peaks used in the interactions(in row 2)
# ------------
PeaksPercentegeInvolved_MANGO_plot_A2=ggplot(MACPET_MACS_Sign_Int_MANGO_data_A2,aes(Merge_Window,100*Peaks_involved_1/Peaks_in,color=factor(Which),linetype=factor(Which)))+
    geom_line(size=0.3)+xlab("MANGO window in bp")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    ggplot2::scale_x_continuous(expand = c(0.005, 0),breaks=seq(0,1000,100))+
    ggplot2::scale_y_continuous(expand = c(0.01,0),breaks=seq(min(round(100*MACPET_MACS_Sign_Int_MANGO_data_A2$Peaks_involved_1/MACPET_MACS_Sign_Int_MANGO_data_A2$Peaks_in)),
                                                              max(round(100*MACPET_MACS_Sign_Int_MANGO_data_A2$Peaks_involved_1/MACPET_MACS_Sign_Int_MANGO_data_A2$Peaks_in)),2))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.15),legend.title=element_blank(),
          legend.text=element_text(size=6),axis.title.y=element_blank())+guides(col = guide_legend(nrow=2))+
    scale_colour_manual(values=c("red","blue","green","brown"))#coloring factors

PeaksPercentegeInvolved_MANGO_plot_A2
save("PeaksPercentegeInvolved_MANGO_plot_A2",file=file.path(DIRsaveA2,"PeaksPercentegeInvolved_MANGO_plot_A2"))


# ------------
# Plot for interactions overlap for 500 window:
# ------------
# find interactions overlap:
# make Granges MACPET:
Anchor1_MACPET_MANGO_A2=GenomicRanges::GRanges(seqnames=MACPETInteractions500_MANGO_A2$chr_from,
                                               ranges=IRanges::IRanges(start=MACPETInteractions500_MANGO_A2$start_from,
                                                                       end=MACPETInteractions500_MANGO_A2$end_from))
Anchor2_MACPET_MANGO_A2=GenomicRanges::GRanges(seqnames=MACPETInteractions500_MANGO_A2$chr_to,
                                               ranges=IRanges::IRanges(start=MACPETInteractions500_MANGO_A2$start_to,
                                                                       end=MACPETInteractions500_MANGO_A2$end_to))
# make interactions
GInteractions_MACPET_MANGO_A2=GInteractions(Anchor1_MACPET_MANGO_A2,Anchor2_MACPET_MANGO_A2)
# make Granges MACS:
Anchor1_MACS_MANGO_A2=GenomicRanges::GRanges(seqnames=MACSInteractions500_MANGO_A2$chr_from,
                                             ranges=IRanges::IRanges(start=MACSInteractions500_MANGO_A2$start_from,
                                                                     end=MACSInteractions500_MANGO_A2$end_from))
Anchor2_MACS_MANGO_A2=GenomicRanges::GRanges(seqnames=MACSInteractions500_MANGO_A2$chr_to,
                                             ranges=IRanges::IRanges(start=MACSInteractions500_MANGO_A2$start_to,
                                                                     end=MACSInteractions500_MANGO_A2$end_to))

GInteractions_MACS_MANGO_A2=GInteractions(Anchor1_MACS_MANGO_A2,Anchor2_MACS_MANGO_A2)


# make calculations for creating the venn data:
library(ggforce)
TotInteractions_MACPET_MANGO_A2=length(GInteractions_MACPET_MANGO_A2)
TotInteractions_MACS_MANGO_A2=length(GInteractions_MACS_MANGO_A2)
TotInteractions_MANGO_A2=TotInteractions_MACPET_MANGO_A2+TotInteractions_MACS_MANGO_A2
TotCommonInteractions_MACStoMACPET_MANGO_A2=length(InteractionSet::findOverlaps(GInteractions_MACPET_MANGO_A2,GInteractions_MACS_MANGO_A2))
# change to 100 scale
TotInteractions_MACPET_MANGO_A2_100=TotInteractions_MACPET_MANGO_A2*100/TotInteractions_MANGO_A2
TotInteractions_MACS_MANGO_A2_100=TotInteractions_MACS_MANGO_A2*100/TotInteractions_MANGO_A2
TotCommonInteractions_MACStoMACPET_MANGO_A2_100=TotCommonInteractions_MACStoMACPET_MANGO_A2*100/TotInteractions_MANGO_A2

VennInteractions_MANGO_A2_data=data.frame(x0=c(0+TotInteractions_MACPET_MANGO_A2_100/2+TotCommonInteractions_MACStoMACPET_MANGO_A2_100/2,#macpet
                                               100-TotInteractions_MACS_MANGO_A2_100/2-TotCommonInteractions_MACStoMACPET_MANGO_A2_100/2,NA),#macs
                                          y0=c(0.5,0.5,NA),
                                          r=c(TotInteractions_MACPET_MANGO_A2_100/2,TotInteractions_MACS_MANGO_A2_100/2,NA),
                                          Totals=c(TotInteractions_MACPET_MANGO_A2,TotInteractions_MACS_MANGO_A2,TotCommonInteractions_MACStoMACPET_MANGO_A2),
                                          color=c("red","blue",NA),which=c("MACPET","MACS",NA))

save("VennInteractions_MANGO_A2_data",file=file.path(DIRsaveA2,"VennInteractions_MANGO_A2_data"))


VennInteractions_Overlap_MANGO_A2_plot=ggplot()+geom_circle(data=VennInteractions_MANGO_A2_data[2,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="blue",alpha=0.6,size=0.2)+
    geom_circle(data=VennInteractions_MANGO_A2_data[1,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="red",alpha=0.6,size=0.2)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),axis.ticks=element_blank(),
          axis.text=element_blank(),axis.title=element_blank())+
    # add numbers for overlaps:
    annotate("text", x=50, y=0.5, label= as.character(VennInteractions_MANGO_A2_data$Totals[3]),size=2.5)+#common
    annotate("text", x=100-VennInteractions_MANGO_A2_data$r[2], y=22, label= as.character(VennInteractions_MANGO_A2_data$Totals[2]-VennInteractions_MANGO_A2_data$Totals[3]),size=2.5)+#macs
    annotate("text", x=0+VennInteractions_MANGO_A2_data$r[1], y=0.5, label= as.character(VennInteractions_MANGO_A2_data$Totals[1]-VennInteractions_MANGO_A2_data$Totals[3]),size=2.5)+#MACPET
    geom_line(data=data.frame(x=rep(100-VennInteractions_MANGO_A2_data$r[2],2),y=c(0,20)),aes(x,y),size=0.3)+
    # add groups:
    annotate("text", x=30, y=30, label= "MACPET",size=3)+#MACPET
    annotate("text", x=70, y=30, label= "MACS",size=3)#MACPET


VennInteractions_Overlap_MANGO_A2_plot
save("VennInteractions_Overlap_MANGO_A2_plot",file=file.path(DIRsaveA2,"VennInteractions_Overlap_MANGO_A2_plot"))

# ------------
# Plot for interactions length for 500 window:
# ------------
scientific_10 <- function(x) {
    y=parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
    if(x[1]==0) y[1]=0

    return(y)
}#for scientific values
# Find interaction Distances:
Interactions_Dist_MACPET_MANGO_A2=pairdist(GInteractions_MACPET_MANGO_A2)
Interactions_Dist_MACS_MANGO_A2=pairdist(GInteractions_MACS_MANGO_A2)
Interaction_Dist_Merged_MANGO_A2=data.frame(Dist=c(Interactions_Dist_MACPET_MANGO_A2,Interactions_Dist_MACS_MANGO_A2),
                                            Which=c(rep("MACPET",length(Interactions_Dist_MACPET_MANGO_A2)),
                                                    rep("MACS",length(Interactions_Dist_MACS_MANGO_A2))))
save("Interaction_Dist_Merged_MANGO_A2",file=file.path(DIRsaveA2,"Interaction_Dist_Merged_MANGO_A2"))


Interaction_Dist_MANGO_Plot_A2=ggplot(Interaction_Dist_Merged_MANGO_A2,aes(Dist,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(label=scientific_10,expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"),axis.title.y=element_blank())+
    xlab("Interactions Distance")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

Interaction_Dist_MANGO_Plot_A2
save("Interaction_Dist_MANGO_Plot_A2",file=file.path(DIRsaveA2,"Interaction_Dist_MANGO_Plot_A2"))

#-------------------------------------------------------
#----------------A3: ENCSR000CAC
#-------------------------------------------------------
# The saving directory
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAC"
Repetition="REP1"
DIRsaveA3=file.path(SaveDir,"ENCSR000CAC/MANGO_interactions")
if(!dir.exists(DIRsaveA3)) dir.create(DIRsaveA3)
DIRdata=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/InteractionsResults_MANGO_MACPETInitials/",sep="")# the data directory
############ Load data MACPET ###############
MACPET_files=list.files(file.path(DIRdata,"MACPET"),pattern="Interactions_")
MACPET_Sign_Int_MANGO_data_A3=data.frame()#for significant interactions
MACPET_HeatMap_mat_MANGO_A3=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACPET_HeatMap_mat_MANGO_A3)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACPET_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACPET_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACPET",folder,"MACPET_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACPET_slopPeak_folder=MACPET_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACPET_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACPET_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACPET_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACPET_slopPeak_folder$chr,
                                                  ranges=IRanges::IRanges(start=MACPET_slopPeak_folder$start,
                                                                          end=MACPET_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACPET",folder,"MACPET.interactions.fdr.mango"))){
        MACPET_Interactions_folder=utils::read.table(file.path(DIRdata,"MACPET",folder,"MACPET.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACPET_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACPET_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACPET_Interactions_folder_update=data.frame(chr=c(MACPET_Interactions_folder$chr_from,
                                                           MACPET_Interactions_folder$chr_to),
                                                     start=c(MACPET_Interactions_folder$start_from,
                                                             MACPET_Interactions_folder$start_to),
                                                     end=c(MACPET_Interactions_folder$end_from,
                                                           MACPET_Interactions_folder$end_to))
        # make granges
        MACPET_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACPET_Interactions_folder_update$chr,
                                                                 ranges=IRanges::IRanges(start=MACPET_Interactions_folder_update$start,
                                                                                         end=MACPET_Interactions_folder_update$end))

        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACPET_slopPeak_folder,MACPET_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACPET_HeatMap_mat_MANGO_A3)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACPET_HeatMap_mat_MANGO_A3[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACPET",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#else no interacions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACPET",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }

    MACPET_Sign_Int_MANGO_data_A3=rbind(MACPET_Sign_Int_MANGO_data_A3,Interactions_folder)
}
# sort by merging window
MACPET_Sign_Int_MANGO_data_A3=MACPET_Sign_Int_MANGO_data_A3[with(MACPET_Sign_Int_MANGO_data_A3,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
MACPET_HeatMap_mat_MANGO_A3=MACPET_HeatMap_mat_MANGO_A3[c(1:(Top_max_used)),]#remove rows
MACPET_HeatMap_mat_MANGO_A3[which(is.na(MACPET_HeatMap_mat_MANGO_A3),arr.ind=TRUE)]=0
save("MACPET_HeatMap_mat_MANGO_A3",file=file.path(DIRsaveA3,"MACPET_HeatMap_mat_MANGO_A3"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
MACPETInteractions500_MANGO_A3=utils::read.table(file.path(DIRdata,"MACPET","Interactions_500","MACPET.interactions.fdr.mango"),stringsAsFactors=FALSE)
colnames(MACPETInteractions500_MANGO_A3)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
save("MACPETInteractions500_MANGO_A3",file=file.path(DIRsaveA3,"MACPETInteractions500_MANGO_A3"))

############ Load data MACS ###############
MACS_files=list.files(file.path(DIRdata,"MACS"),pattern="Interactions_")
MACS_0.05_Sign_Int_MANGO_data_A3=data.frame()
MACS_HeatMap_mat_MANGO_A3=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A3)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A3)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A3[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }

        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.05",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.05",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }

    MACS_0.05_Sign_Int_MANGO_data_A3=rbind(MACS_0.05_Sign_Int_MANGO_data_A3,Interactions_folder)
}
# sort by merging window
MACS_0.05_Sign_Int_MANGO_data_A3=MACS_0.05_Sign_Int_MANGO_data_A3[with(MACS_0.05_Sign_Int_MANGO_data_A3,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
MACS_HeatMap_mat_MANGO_A3=MACS_HeatMap_mat_MANGO_A3[c(1:(Top_max_used)),]#remove rows
MACS_HeatMap_mat_MANGO_A3[which(is.na(MACS_HeatMap_mat_MANGO_A3),arr.ind=TRUE)]=0
save("MACS_HeatMap_mat_MANGO_A3",file=file.path(DIRsaveA3,"MACS_HeatMap_mat_MANGO_A3"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
MACSInteractions500_MANGO_A3=utils::read.table(file.path(DIRdata,"MACS","Interactions_500","MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
colnames(MACSInteractions500_MANGO_A3)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
save("MACSInteractions500_MANGO_A3",file=file.path(DIRsaveA3,"MACSInteractions500_MANGO_A3"))

############ Load data MACS 0.01 ###############
MACS_0.01_files=list.files(file.path(DIRdata,"MACS_0.01"),pattern="Interactions_")
MACS_0.01_Sign_Int_MANGO_data_A3=data.frame()
MACS_HeatMap_mat_MANGO_A3=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A3)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_0.01_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS_0.01",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS_0.01",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS_0.01",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A3)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A3[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.01",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.01",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }
    
    MACS_0.01_Sign_Int_MANGO_data_A3=rbind(MACS_0.01_Sign_Int_MANGO_data_A3,Interactions_folder)
}
# sort by merging window
MACS_0.01_Sign_Int_MANGO_data_A3=MACS_0.01_Sign_Int_MANGO_data_A3[with(MACS_0.01_Sign_Int_MANGO_data_A3,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
# MACS_HeatMap_mat_MANGO_A3=MACS_HeatMap_mat_MANGO_A3[c(1:(Top_max_used)),]#remove rows
# MACS_HeatMap_mat_MANGO_A3[which(is.na(MACS_HeatMap_mat_MANGO_A3),arr.ind=TRUE)]=0
# save("MACS_HeatMap_mat_MANGO_A3",file=file.path(DIRsaveA3,"MACS_HeatMap_mat_MANGO_A3"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
# MACSInteractions500_MANGO_A3=utils::read.table(file.path(DIRdata,"MACS","Interactions_500","MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
# colnames(MACSInteractions500_MANGO_A3)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
# save("MACSInteractions500_MANGO_A3",file=file.path(DIRsaveA3,"MACSInteractions500_MANGO_A3"))
############ Load data MACS ###############
MACS_N_files=list.files(file.path(DIRdata,"MACS_N"),pattern="Interactions_")
MACS_N_Sign_Int_MANGO_data_A3=data.frame()
MACS_HeatMap_mat_MANGO_A3=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A3)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_N_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS_N",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS_N",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS_N",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A3)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A3[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS N",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS N",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }
    
    MACS_N_Sign_Int_MANGO_data_A3=rbind(MACS_N_Sign_Int_MANGO_data_A3,Interactions_folder)
}
# sort by merging window
MACS_N_Sign_Int_MANGO_data_A3=MACS_N_Sign_Int_MANGO_data_A3[with(MACS_N_Sign_Int_MANGO_data_A3,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
# MACS_HeatMap_mat_MANGO_A3=MACS_HeatMap_mat_MANGO_A3[c(1:(Top_max_used)),]#remove rows
# MACS_HeatMap_mat_MANGO_A3[which(is.na(MACS_HeatMap_mat_MANGO_A3),arr.ind=TRUE)]=0
# save("MACS_HeatMap_mat_MANGO_A3",file=file.path(DIRsaveA3,"MACS_HeatMap_mat_MANGO_A3"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
# MACSInteractions500_MANGO_A3=utils::read.table(file.path(DIRdata,"MACS","Interactions_500","MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
# colnames(MACSInteractions500_MANGO_A3)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
# save("MACSInteractions500_MANGO_A3",file=file.path(DIRsaveA3,"MACSInteractions500_MANGO_A3"))


############ Merge and make plot ###############
MACPET_MACS_Sign_Int_MANGO_data_A3=rbind(MACPET_Sign_Int_MANGO_data_A3,MACS_0.05_Sign_Int_MANGO_data_A3,MACS_0.01_Sign_Int_MANGO_data_A3,MACS_N_Sign_Int_MANGO_data_A3)
save("MACPET_MACS_Sign_Int_MANGO_data_A3",file=file.path(DIRsaveA3,"MACPET_MACS_Sign_Int_MANGO_data_A3"))

# ------------
# significant interactions plot (in row 1):
# ------------
SigInteractions_MANGO_plot_A3=ggplot(MACPET_MACS_Sign_Int_MANGO_data_A3,aes(Merge_Window,Total_significant_int,color=factor(Which),linetype=factor(Which)))+
    geom_line(size=0.3)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.x=element_blank(),axis.title.y=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.005, 0),breaks=seq(0,1000,100))+
    ggplot2::scale_y_continuous(expand = c(0.01,0),breaks=seq(min(MACPET_MACS_Sign_Int_MANGO_data_A3$Total_significant_int),max(MACPET_MACS_Sign_Int_MANGO_data_A3$Total_significant_int),500))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.15),legend.title=element_blank(),
          legend.text=element_text(size=6))+guides(col = guide_legend(nrow=2))+
    scale_colour_manual(values=c("red","blue","green","brown"))#coloring factors

SigInteractions_MANGO_plot_A3
save("SigInteractions_MANGO_plot_A3",file=file.path(DIRsaveA3,"SigInteractions_MANGO_plot_A3"))


# ------------
# plots for percentage of peaks inn and peaks used in the interactions(in row 2)
# ------------
PeaksPercentegeInvolved_MANGO_plot_A3=ggplot(MACPET_MACS_Sign_Int_MANGO_data_A3,aes(Merge_Window,100*Peaks_involved_1/Peaks_in,color=factor(Which),linetype=factor(Which)))+
    geom_line(size=0.3)+xlab(" ")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"),axis.title.y=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.005, 0),breaks=seq(0,1000,100))+
    ggplot2::scale_y_continuous(expand = c(0.01,0),breaks=seq(min(round(100*MACPET_MACS_Sign_Int_MANGO_data_A3$Peaks_involved_1/MACPET_MACS_Sign_Int_MANGO_data_A3$Peaks_in)),
                                                              max(round(100*MACPET_MACS_Sign_Int_MANGO_data_A3$Peaks_involved_1/MACPET_MACS_Sign_Int_MANGO_data_A3$Peaks_in)),2))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.25),legend.title=element_blank(),
          legend.text=element_text(size=6))+guides(col = guide_legend(nrow=2))+
    scale_colour_manual(values=c("red","blue","green","brown"))#coloring factors

PeaksPercentegeInvolved_MANGO_plot_A3
save("PeaksPercentegeInvolved_MANGO_plot_A3",file=file.path(DIRsaveA3,"PeaksPercentegeInvolved_MANGO_plot_A3"))


# ------------
# Plot for interactions overlap for 500 window:
# ------------
# find interactions overlap:
# make Granges MACPET:
Anchor1_MACPET_MANGO_A3=GenomicRanges::GRanges(seqnames=MACPETInteractions500_MANGO_A3$chr_from,
                                               ranges=IRanges::IRanges(start=MACPETInteractions500_MANGO_A3$start_from,
                                                                       end=MACPETInteractions500_MANGO_A3$end_from))
Anchor2_MACPET_MANGO_A3=GenomicRanges::GRanges(seqnames=MACPETInteractions500_MANGO_A3$chr_to,
                                               ranges=IRanges::IRanges(start=MACPETInteractions500_MANGO_A3$start_to,
                                                                       end=MACPETInteractions500_MANGO_A3$end_to))
# make interactions
GInteractions_MACPET_MANGO_A3=GInteractions(Anchor1_MACPET_MANGO_A3,Anchor2_MACPET_MANGO_A3)
# make Granges MACS:
Anchor1_MACS_MANGO_A3=GenomicRanges::GRanges(seqnames=MACSInteractions500_MANGO_A3$chr_from,
                                             ranges=IRanges::IRanges(start=MACSInteractions500_MANGO_A3$start_from,
                                                                     end=MACSInteractions500_MANGO_A3$end_from))
Anchor2_MACS_MANGO_A3=GenomicRanges::GRanges(seqnames=MACSInteractions500_MANGO_A3$chr_to,
                                             ranges=IRanges::IRanges(start=MACSInteractions500_MANGO_A3$start_to,
                                                                     end=MACSInteractions500_MANGO_A3$end_to))

GInteractions_MACS_MANGO_A3=GInteractions(Anchor1_MACS_MANGO_A3,Anchor2_MACS_MANGO_A3)


# make calculations for creating the venn data:
library(ggforce)
TotInteractions_MACPET_MANGO_A3=length(GInteractions_MACPET_MANGO_A3)
TotInteractions_MACS_MANGO_A3=length(GInteractions_MACS_MANGO_A3)
TotInteractions_MANGO_A3=TotInteractions_MACPET_MANGO_A3+TotInteractions_MACS_MANGO_A3
TotCommonInteractions_MACStoMACPET_MANGO_A3=length(InteractionSet::findOverlaps(GInteractions_MACPET_MANGO_A3,GInteractions_MACS_MANGO_A3))
# change to 100 scale
TotInteractions_MACPET_MANGO_A3_100=TotInteractions_MACPET_MANGO_A3*100/TotInteractions_MANGO_A3
TotInteractions_MACS_MANGO_A3_100=TotInteractions_MACS_MANGO_A3*100/TotInteractions_MANGO_A3
TotCommonInteractions_MACStoMACPET_MANGO_A3_100=TotCommonInteractions_MACStoMACPET_MANGO_A3*100/TotInteractions_MANGO_A3

VennInteractions_MANGO_A3_data=data.frame(x0=c(0+TotInteractions_MACPET_MANGO_A3_100/2+TotCommonInteractions_MACStoMACPET_MANGO_A3_100/2,#macpet
                                               100-TotInteractions_MACS_MANGO_A3_100/2-TotCommonInteractions_MACStoMACPET_MANGO_A3_100/2,NA),#macs
                                          y0=c(0.5,0.5,NA),
                                          r=c(TotInteractions_MACPET_MANGO_A3_100/2,TotInteractions_MACS_MANGO_A3_100/2,NA),
                                          Totals=c(TotInteractions_MACPET_MANGO_A3,TotInteractions_MACS_MANGO_A3,TotCommonInteractions_MACStoMACPET_MANGO_A3),
                                          color=c("red","blue",NA),which=c("MACPET","MACS",NA))

save("VennInteractions_MANGO_A3_data",file=file.path(DIRsaveA3,"VennInteractions_MANGO_A3_data"))


VennInteractions_Overlap_MANGO_A3_plot=ggplot()+geom_circle(data=VennInteractions_MANGO_A3_data[2,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="blue",alpha=0.6,size=0.2)+
    geom_circle(data=VennInteractions_MANGO_A3_data[1,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="red",alpha=0.6,size=0.2)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),axis.ticks=element_blank(),
          axis.text=element_blank(),axis.title=element_blank())+
    # add numbers for overlaps:
    annotate("text", x=55, y=0.5, label= as.character(VennInteractions_MANGO_A3_data$Totals[3]),size=2.5)+#common
    annotate("text", x=100-VennInteractions_MANGO_A3_data$r[2], y=22, label= as.character(VennInteractions_MANGO_A3_data$Totals[2]-VennInteractions_MANGO_A3_data$Totals[3]),size=2.5)+#macs
    annotate("text", x=0+VennInteractions_MANGO_A3_data$r[1], y=0.5, label= as.character(VennInteractions_MANGO_A3_data$Totals[1]-VennInteractions_MANGO_A3_data$Totals[3]),size=2.5)+#MACPET
    geom_line(data=data.frame(x=rep(100-VennInteractions_MANGO_A3_data$r[2],2),y=c(0,20)),aes(x,y),size=0.3)+
    # add groups:
    annotate("text", x=30, y=30, label= "MACPET",size=3)+#MACPET
    annotate("text", x=70, y=30, label= "MACS",size=3)#MACPET


VennInteractions_Overlap_MANGO_A3_plot
save("VennInteractions_Overlap_MANGO_A3_plot",file=file.path(DIRsaveA3,"VennInteractions_Overlap_MANGO_A3_plot"))

# ------------
# Plot for interactions length for 500 window:
# ------------
scientific_10 <- function(x) {
    y=parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
    if(x[1]==0) y[1]=0

    return(y)
}#for scientific values
# Find interaction Distances:
Interactions_Dist_MACPET_MANGO_A3=pairdist(GInteractions_MACPET_MANGO_A3)
Interactions_Dist_MACS_MANGO_A3=pairdist(GInteractions_MACS_MANGO_A3)
Interaction_Dist_Merged_MANGO_A3=data.frame(Dist=c(Interactions_Dist_MACPET_MANGO_A3,Interactions_Dist_MACS_MANGO_A3),
                                            Which=c(rep("MACPET",length(Interactions_Dist_MACPET_MANGO_A3)),
                                                    rep("MACS",length(Interactions_Dist_MACS_MANGO_A3))))
save("Interaction_Dist_Merged_MANGO_A3",file=file.path(DIRsaveA3,"Interaction_Dist_Merged_MANGO_A3"))


Interaction_Dist_MANGO_Plot_A3=ggplot(Interaction_Dist_Merged_MANGO_A3,aes(Dist,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(label=scientific_10,expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0),breaks=c(250000,500000,750000,950000))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"),axis.title.y=element_blank())+
    xlab("")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

Interaction_Dist_MANGO_Plot_A3
save("Interaction_Dist_MANGO_Plot_A3",file=file.path(DIRsaveA3,"Interaction_Dist_MANGO_Plot_A3"))



#-------------------------------------------------------
#----------------A4: ENCSR000FDD
#-------------------------------------------------------
# The saving directory
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000FDD"
Repetition="REP1"
DIRsaveA4=file.path(SaveDir,"ENCSR000FDD/MANGO_interactions")
if(!dir.exists(DIRsaveA4)) dir.create(DIRsaveA4)
DIRdata=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/InteractionsResults_MANGO_MACPETInitials/",sep="")# the data directory
############ Load data MACPET ###############
MACPET_files=list.files(file.path(DIRdata,"MACPET"),pattern="Interactions_")
MACPET_Sign_Int_MANGO_data_A4=data.frame()#for significant interactions
MACPET_HeatMap_mat_MANGO_A4=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACPET_HeatMap_mat_MANGO_A4)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACPET_files){#CHANGES
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACPET_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACPET",folder,"MACPET_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACPET_slopPeak_folder=MACPET_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACPET_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACPET_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACPET_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACPET_slopPeak_folder$chr,
                                                  ranges=IRanges::IRanges(start=MACPET_slopPeak_folder$start,
                                                                          end=MACPET_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACPET",folder,"MACPET.interactions.fdr.mango"))){
        MACPET_Interactions_folder=utils::read.table(file.path(DIRdata,"MACPET",folder,"MACPET.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACPET_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACPET_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACPET_Interactions_folder_update=data.frame(chr=c(MACPET_Interactions_folder$chr_from,
                                                           MACPET_Interactions_folder$chr_to),
                                                     start=c(MACPET_Interactions_folder$start_from,
                                                             MACPET_Interactions_folder$start_to),
                                                     end=c(MACPET_Interactions_folder$end_from,
                                                           MACPET_Interactions_folder$end_to))
        # make granges
        MACPET_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACPET_Interactions_folder_update$chr,
                                                                 ranges=IRanges::IRanges(start=MACPET_Interactions_folder_update$start,
                                                                                         end=MACPET_Interactions_folder_update$end))

        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACPET_slopPeak_folder,MACPET_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACPET_HeatMap_mat_MANGO_A4)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACPET_HeatMap_mat_MANGO_A4[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACPET",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#else no interacions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACPET",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }

    MACPET_Sign_Int_MANGO_data_A4=rbind(MACPET_Sign_Int_MANGO_data_A4,Interactions_folder)
}
# sort by merging window
MACPET_Sign_Int_MANGO_data_A4=MACPET_Sign_Int_MANGO_data_A4[with(MACPET_Sign_Int_MANGO_data_A4,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
MACPET_HeatMap_mat_MANGO_A4=MACPET_HeatMap_mat_MANGO_A4[c(1:(Top_max_used)),]#remove rows
MACPET_HeatMap_mat_MANGO_A4[which(is.na(MACPET_HeatMap_mat_MANGO_A4),arr.ind=TRUE)]=0
save("MACPET_HeatMap_mat_MANGO_A4",file=file.path(DIRsaveA4,"MACPET_HeatMap_mat_MANGO_A4"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
MACPETInteractions500_MANGO_A4=utils::read.table(file.path(DIRdata,"MACPET","Interactions_500","MACPET.interactions.fdr.mango"),stringsAsFactors=FALSE)
colnames(MACPETInteractions500_MANGO_A4)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
save("MACPETInteractions500_MANGO_A4",file=file.path(DIRsaveA4,"MACPETInteractions500_MANGO_A4"))

############ Load data MACS ###############
MACS_files=list.files(file.path(DIRdata,"MACS"),pattern="Interactions_")
MACS_0.05_Sign_Int_MANGO_data_A4=data.frame()
MACS_HeatMap_mat_MANGO_A4=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A4)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A4)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A4[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }

        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.05",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.05",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }

    MACS_0.05_Sign_Int_MANGO_data_A4=rbind(MACS_0.05_Sign_Int_MANGO_data_A4,Interactions_folder)
}
# sort by merging window
MACS_0.05_Sign_Int_MANGO_data_A4=MACS_0.05_Sign_Int_MANGO_data_A4[with(MACS_0.05_Sign_Int_MANGO_data_A4,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
MACS_HeatMap_mat_MANGO_A4=MACS_HeatMap_mat_MANGO_A4[c(1:(Top_max_used)),]#remove rows
MACS_HeatMap_mat_MANGO_A4[which(is.na(MACS_HeatMap_mat_MANGO_A4),arr.ind=TRUE)]=0
save("MACS_HeatMap_mat_MANGO_A4",file=file.path(DIRsaveA4,"MACS_HeatMap_mat_MANGO_A4"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
MACSInteractions500_MANGO_A4=utils::read.table(file.path(DIRdata,"MACS","Interactions_500","MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
colnames(MACSInteractions500_MANGO_A4)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
save("MACSInteractions500_MANGO_A4",file=file.path(DIRsaveA4,"MACSInteractions500_MANGO_A4"))

############ Load data MACS 0.01 ###############
MACS_0.01_files=list.files(file.path(DIRdata,"MACS_0.01"),pattern="Interactions_")
MACS_0.01_Sign_Int_MANGO_data_A4=data.frame()
MACS_HeatMap_mat_MANGO_A4=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A4)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_0.01_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS_0.01",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS_0.01",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS_0.01",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A4)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A4[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.01",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.01",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }
    
    MACS_0.01_Sign_Int_MANGO_data_A4=rbind(MACS_0.01_Sign_Int_MANGO_data_A4,Interactions_folder)
}
# sort by merging window
MACS_0.01_Sign_Int_MANGO_data_A4=MACS_0.01_Sign_Int_MANGO_data_A4[with(MACS_0.01_Sign_Int_MANGO_data_A4,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
# MACS_HeatMap_mat_MANGO_A4=MACS_HeatMap_mat_MANGO_A4[c(1:(Top_max_used)),]#remove rows
# MACS_HeatMap_mat_MANGO_A4[which(is.na(MACS_HeatMap_mat_MANGO_A4),arr.ind=TRUE)]=0
# save("MACS_HeatMap_mat_MANGO_A4",file=file.path(DIRsaveA4,"MACS_HeatMap_mat_MANGO_A4"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
# MACSInteractions500_MANGO_A4=utils::read.table(file.path(DIRdata,"MACS","Interactions_500","MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
# colnames(MACSInteractions500_MANGO_A4)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
# save("MACSInteractions500_MANGO_A4",file=file.path(DIRsaveA4,"MACSInteractions500_MANGO_A4"))

############ Load data MACS N ###############
MACS_N_files=list.files(file.path(DIRdata,"MACS_N"),pattern="Interactions_")
MACS_N_Sign_Int_MANGO_data_A4=data.frame()
MACS_HeatMap_mat_MANGO_A4=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A4)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_N_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS_N",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS_N",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS_N",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A4)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A4[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS N",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS N",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }
    
    MACS_N_Sign_Int_MANGO_data_A4=rbind(MACS_N_Sign_Int_MANGO_data_A4,Interactions_folder)
}
# sort by merging window
MACS_N_Sign_Int_MANGO_data_A4=MACS_N_Sign_Int_MANGO_data_A4[with(MACS_N_Sign_Int_MANGO_data_A4,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
# MACS_HeatMap_mat_MANGO_A4=MACS_HeatMap_mat_MANGO_A4[c(1:(Top_max_used)),]#remove rows
# MACS_HeatMap_mat_MANGO_A4[which(is.na(MACS_HeatMap_mat_MANGO_A4),arr.ind=TRUE)]=0
# save("MACS_HeatMap_mat_MANGO_A4",file=file.path(DIRsaveA4,"MACS_HeatMap_mat_MANGO_A4"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
# MACSInteractions500_MANGO_A4=utils::read.table(file.path(DIRdata,"MACS","Interactions_500","MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
# colnames(MACSInteractions500_MANGO_A4)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
# save("MACSInteractions500_MANGO_A4",file=file.path(DIRsaveA4,"MACSInteractions500_MANGO_A4"))


############ Merge and make plot ###############
MACPET_MACS_Sign_Int_MANGO_data_A4=rbind(MACPET_Sign_Int_MANGO_data_A4,MACS_0.05_Sign_Int_MANGO_data_A4,MACS_0.01_Sign_Int_MANGO_data_A4,MACS_N_Sign_Int_MANGO_data_A4)
save("MACPET_MACS_Sign_Int_MANGO_data_A4",file=file.path(DIRsaveA4,"MACPET_MACS_Sign_Int_MANGO_data_A4"))

# ------------
# significant interactions plot (in row 1):
# ------------
SigInteractions_MANGO_plot_A4=ggplot(MACPET_MACS_Sign_Int_MANGO_data_A4,aes(Merge_Window,Total_significant_int,color=factor(Which),linetype=factor(Which)))+
    geom_line(size=0.3)+ylab("Significant Interactions")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.x=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.005, 0),breaks=seq(0,1000,100))+
    ggplot2::scale_y_continuous(expand = c(0.01,0),breaks=seq(min(MACPET_MACS_Sign_Int_MANGO_data_A4$Total_significant_int),
                                                              max(MACPET_MACS_Sign_Int_MANGO_data_A4$Total_significant_int),10))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(1.07,.4),legend.title=element_blank(),
          legend.text=element_text(size=6),axis.title.x=element_blank())+guides(col = guide_legend(nrow=2))+
    scale_colour_manual(values=c("red","blue","green","brown"))#coloring factors

SigInteractions_MANGO_plot_A4
save("SigInteractions_MANGO_plot_A4",file=file.path(DIRsaveA4,"SigInteractions_MANGO_plot_A4"))


# ------------
# plots for percentage of peaks inn and peaks used in the interactions(in row 2)
# ------------
PeaksPercentegeInvolved_MANGO_plot_A4=ggplot(MACPET_MACS_Sign_Int_MANGO_data_A4,aes(Merge_Window,100*Peaks_involved_1/Peaks_in,color=factor(Which),linetype=factor(Which)))+
    geom_line(size=0.3)+ylab("Peaks used (%)")+xlab(" ")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    ggplot2::scale_x_continuous(expand = c(0.005, 0),breaks=seq(0,1000,100))+
    ggplot2::scale_y_continuous(expand = c(0.01,0),breaks=seq(min(round(100*MACPET_MACS_Sign_Int_MANGO_data_A4$Peaks_involved_1/MACPET_MACS_Sign_Int_MANGO_data_A4$Peaks_in)),
                                                              max(round(100*MACPET_MACS_Sign_Int_MANGO_data_A4$Peaks_involved_1/MACPET_MACS_Sign_Int_MANGO_data_A4$Peaks_in)),0.1))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.20),legend.title=element_blank(),
          legend.text=element_text(size=6))+guides(col = guide_legend(nrow=2))+
    scale_colour_manual(values=c("red","blue","green","brown"))#coloring factors

PeaksPercentegeInvolved_MANGO_plot_A4
save("PeaksPercentegeInvolved_MANGO_plot_A4",file=file.path(DIRsaveA4,"PeaksPercentegeInvolved_MANGO_plot_A4"))


# ------------
# Plot for interactions overlap for 500 window:
# ------------
# find interactions overlap:
# make Granges MACPET:
Anchor1_MACPET_MANGO_A4=GenomicRanges::GRanges(seqnames=MACPETInteractions500_MANGO_A4$chr_from,
                                               ranges=IRanges::IRanges(start=MACPETInteractions500_MANGO_A4$start_from,
                                                                       end=MACPETInteractions500_MANGO_A4$end_from))
Anchor2_MACPET_MANGO_A4=GenomicRanges::GRanges(seqnames=MACPETInteractions500_MANGO_A4$chr_to,
                                               ranges=IRanges::IRanges(start=MACPETInteractions500_MANGO_A4$start_to,
                                                                       end=MACPETInteractions500_MANGO_A4$end_to))
# make interactions
GInteractions_MACPET_MANGO_A4=GInteractions(Anchor1_MACPET_MANGO_A4,Anchor2_MACPET_MANGO_A4)
# make Granges MACS:
Anchor1_MACS_MANGO_A4=GenomicRanges::GRanges(seqnames=MACSInteractions500_MANGO_A4$chr_from,
                                             ranges=IRanges::IRanges(start=MACSInteractions500_MANGO_A4$start_from,
                                                                     end=MACSInteractions500_MANGO_A4$end_from))
Anchor2_MACS_MANGO_A4=GenomicRanges::GRanges(seqnames=MACSInteractions500_MANGO_A4$chr_to,
                                             ranges=IRanges::IRanges(start=MACSInteractions500_MANGO_A4$start_to,
                                                                     end=MACSInteractions500_MANGO_A4$end_to))

GInteractions_MACS_MANGO_A4=GInteractions(Anchor1_MACS_MANGO_A4,Anchor2_MACS_MANGO_A4)


# make calculations for creating the venn data:
library(ggforce)
TotInteractions_MACPET_MANGO_A4=length(GInteractions_MACPET_MANGO_A4)
TotInteractions_MACS_MANGO_A4=length(GInteractions_MACS_MANGO_A4)
TotInteractions_MANGO_A4=TotInteractions_MACPET_MANGO_A4+TotInteractions_MACS_MANGO_A4
TotCommonInteractions_MACStoMACPET_MANGO_A4=length(InteractionSet::findOverlaps(GInteractions_MACPET_MANGO_A4,GInteractions_MACS_MANGO_A4))
# change to 100 scale
TotInteractions_MACPET_MANGO_A4_100=TotInteractions_MACPET_MANGO_A4*100/TotInteractions_MANGO_A4
TotInteractions_MACS_MANGO_A4_100=TotInteractions_MACS_MANGO_A4*100/TotInteractions_MANGO_A4
TotCommonInteractions_MACStoMACPET_MANGO_A4_100=TotCommonInteractions_MACStoMACPET_MANGO_A4*100/TotInteractions_MANGO_A4

VennInteractions_MANGO_A4_data=data.frame(x0=c(0+TotInteractions_MACPET_MANGO_A4_100/2+TotCommonInteractions_MACStoMACPET_MANGO_A4_100/2,#macpet
                                               100-TotInteractions_MACS_MANGO_A4_100/2-TotCommonInteractions_MACStoMACPET_MANGO_A4_100/2,NA),#macs
                                          y0=c(0.5,0.5,NA),
                                          r=c(TotInteractions_MACPET_MANGO_A4_100/2,TotInteractions_MACS_MANGO_A4_100/2,NA),
                                          Totals=c(TotInteractions_MACPET_MANGO_A4,TotInteractions_MACS_MANGO_A4,TotCommonInteractions_MACStoMACPET_MANGO_A4),
                                          color=c("red","blue",NA),which=c("MACPET","MACS",NA))

save("VennInteractions_MANGO_A4_data",file=file.path(DIRsaveA4,"VennInteractions_MANGO_A4_data"))


VennInteractions_Overlap_MANGO_A4_plot=ggplot()+geom_circle(data=VennInteractions_MANGO_A4_data[2,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="blue",alpha=0.6,size=0.2)+
    geom_circle(data=VennInteractions_MANGO_A4_data[1,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="red",alpha=0.6,size=0.2)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),axis.ticks=element_blank(),
          axis.text=element_blank(),axis.title=element_blank())+
    # add numbers for overlaps:
    annotate("text", x=63, y=0.5, label= as.character(VennInteractions_MANGO_A4_data$Totals[3]),size=2.5)+#common
    annotate("text", x=100-VennInteractions_MANGO_A4_data$r[2], y=0.5, label= as.character(VennInteractions_MANGO_A4_data$Totals[2]-VennInteractions_MANGO_A4_data$Totals[3]),size=2.5)+#macs
    annotate("text", x=0+VennInteractions_MANGO_A4_data$r[1], y=0.5, label= as.character(VennInteractions_MANGO_A4_data$Totals[1]-VennInteractions_MANGO_A4_data$Totals[3]),size=2.5)+#MACPET
    # add groups:
    annotate("text", x=20, y=40, label= "MACPET",size=3)+#MACPET
    annotate("text", x=75, y=40, label= "MACS",size=3)#MACPET


VennInteractions_Overlap_MANGO_A4_plot
save("VennInteractions_Overlap_MANGO_A4_plot",file=file.path(DIRsaveA4,"VennInteractions_Overlap_MANGO_A4_plot"))

# ------------
# Plot for interactions length for 500 window:
# ------------
scientific_10 <- function(x) {
    y=parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
    if(x[1]==0) y[1]=0

    return(y)
}#for scientific values
# Find interaction Distances:
Interactions_Dist_MACPET_MANGO_A4=pairdist(GInteractions_MACPET_MANGO_A4)
Interactions_Dist_MACS_MANGO_A4=pairdist(GInteractions_MACS_MANGO_A4)
Interaction_Dist_Merged_MANGO_A4=data.frame(Dist=c(Interactions_Dist_MACPET_MANGO_A4,Interactions_Dist_MACS_MANGO_A4),
                                            Which=c(rep("MACPET",length(Interactions_Dist_MACPET_MANGO_A4)),
                                                    rep("MACS",length(Interactions_Dist_MACS_MANGO_A4))))
save("Interaction_Dist_Merged_MANGO_A4",file=file.path(DIRsaveA4,"Interaction_Dist_Merged_MANGO_A4"))


Interaction_Dist_MANGO_Plot_A4=ggplot(Interaction_Dist_Merged_MANGO_A4,aes(Dist,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(label=scientific_10,expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    xlab("")+
    ylab("Frequency")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

Interaction_Dist_MANGO_Plot_A4
save("Interaction_Dist_MANGO_Plot_A4",file=file.path(DIRsaveA4,"Interaction_Dist_MANGO_Plot_A4"))

#-------------------------------------------------------
#----------------A5: ENCSR000FDG or Not_specified_1
#-------------------------------------------------------
# The saving directory
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000FDG"
Repetition="REP1"
DIRsaveA5=file.path(SaveDir,"ENCSR000FDG/MANGO_interactions")
if(!dir.exists(DIRsaveA5)) dir.create(DIRsaveA5)
DIRdata=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/InteractionsResults_MANGO_MACPETInitials/",sep="")# the data directory
############ Load data MACPET ###############
MACPET_files=list.files(file.path(DIRdata,"MACPET"),pattern="Interactions_")
MACPET_Sign_Int_MANGO_data_A5=data.frame()#for significant interactions
MACPET_HeatMap_mat_MANGO_A5=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACPET_HeatMap_mat_MANGO_A5)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACPET_files){#CHANGES
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACPET_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACPET",folder,"MACPET_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACPET_slopPeak_folder=MACPET_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACPET_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACPET_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACPET_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACPET_slopPeak_folder$chr,
                                                  ranges=IRanges::IRanges(start=MACPET_slopPeak_folder$start,
                                                                          end=MACPET_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACPET",folder,"MACPET.interactions.fdr.mango"))){
        MACPET_Interactions_folder=utils::read.table(file.path(DIRdata,"MACPET",folder,"MACPET.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACPET_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACPET_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACPET_Interactions_folder_update=data.frame(chr=c(MACPET_Interactions_folder$chr_from,
                                                           MACPET_Interactions_folder$chr_to),
                                                     start=c(MACPET_Interactions_folder$start_from,
                                                             MACPET_Interactions_folder$start_to),
                                                     end=c(MACPET_Interactions_folder$end_from,
                                                           MACPET_Interactions_folder$end_to))
        # make granges
        MACPET_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACPET_Interactions_folder_update$chr,
                                                                 ranges=IRanges::IRanges(start=MACPET_Interactions_folder_update$start,
                                                                                         end=MACPET_Interactions_folder_update$end))

        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACPET_slopPeak_folder,MACPET_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACPET_HeatMap_mat_MANGO_A5)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACPET_HeatMap_mat_MANGO_A5[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACPET",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#else no interacions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACPET",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }

    MACPET_Sign_Int_MANGO_data_A5=rbind(MACPET_Sign_Int_MANGO_data_A5,Interactions_folder)
}
# sort by merging window
MACPET_Sign_Int_MANGO_data_A5=MACPET_Sign_Int_MANGO_data_A5[with(MACPET_Sign_Int_MANGO_data_A5,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
MACPET_HeatMap_mat_MANGO_A5=MACPET_HeatMap_mat_MANGO_A5[c(1:(Top_max_used)),]#remove rows
MACPET_HeatMap_mat_MANGO_A5[which(is.na(MACPET_HeatMap_mat_MANGO_A5),arr.ind=TRUE)]=0
save("MACPET_HeatMap_mat_MANGO_A5",file=file.path(DIRsaveA5,"MACPET_HeatMap_mat_MANGO_A5"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
MACPETInteractions500_MANGO_A5=utils::read.table(file.path(DIRdata,"MACPET","Interactions_500","MACPET.interactions.fdr.mango"),stringsAsFactors=FALSE)
colnames(MACPETInteractions500_MANGO_A5)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
save("MACPETInteractions500_MANGO_A5",file=file.path(DIRsaveA5,"MACPETInteractions500_MANGO_A5"))

############ Load data MACS ###############
MACS_files=list.files(file.path(DIRdata,"MACS"),pattern="Interactions_")
MACS_0.05_Sign_Int_MANGO_data_A5=data.frame()
MACS_HeatMap_mat_MANGO_A5=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A5)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A5)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A5[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }

        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.05",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.05",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }

    MACS_0.05_Sign_Int_MANGO_data_A5=rbind(MACS_0.05_Sign_Int_MANGO_data_A5,Interactions_folder)
}
# sort by merging window
MACS_0.05_Sign_Int_MANGO_data_A5=MACS_0.05_Sign_Int_MANGO_data_A5[with(MACS_0.05_Sign_Int_MANGO_data_A5,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
MACS_HeatMap_mat_MANGO_A5=MACS_HeatMap_mat_MANGO_A5[c(1:(Top_max_used)),]#remove rows
MACS_HeatMap_mat_MANGO_A5[which(is.na(MACS_HeatMap_mat_MANGO_A5),arr.ind=TRUE)]=0
save("MACS_HeatMap_mat_MANGO_A5",file=file.path(DIRsaveA5,"MACS_HeatMap_mat_MANGO_A5"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
MACSInteractions500_MANGO_A5=utils::read.table(file.path(DIRdata,"MACS","Interactions_500","MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
colnames(MACSInteractions500_MANGO_A5)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
save("MACSInteractions500_MANGO_A5",file=file.path(DIRsaveA5,"MACSInteractions500_MANGO_A5"))

############ Load data MACS 0.01 ###############
MACS_0.01_files=list.files(file.path(DIRdata,"MACS_0.01"),pattern="Interactions_")
MACS_0.01_Sign_Int_MANGO_data_A5=data.frame()
MACS_HeatMap_mat_MANGO_A5=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A5)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_0.01_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS_0.01",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS_0.01",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS_0.01",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A5)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A5[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.01",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.01",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }
    
    MACS_0.01_Sign_Int_MANGO_data_A5=rbind(MACS_0.01_Sign_Int_MANGO_data_A5,Interactions_folder)
}
# sort by merging window
MACS_0.01_Sign_Int_MANGO_data_A5=MACS_0.01_Sign_Int_MANGO_data_A5[with(MACS_0.01_Sign_Int_MANGO_data_A5,order(Merge_Window)),]#used for significant interactions plot

############ Load data MACS N ###############
MACS_N_files=list.files(file.path(DIRdata,"MACS_N"),pattern="Interactions_")
MACS_N_Sign_Int_MANGO_data_A5=data.frame()
MACS_HeatMap_mat_MANGO_A5=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A5)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_N_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS_N",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS_N",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS_N",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A5)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A5[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS N",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS N",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }
    
    MACS_N_Sign_Int_MANGO_data_A5=rbind(MACS_N_Sign_Int_MANGO_data_A5,Interactions_folder)
}
# sort by merging window
MACS_N_Sign_Int_MANGO_data_A5=MACS_N_Sign_Int_MANGO_data_A5[with(MACS_N_Sign_Int_MANGO_data_A5,order(Merge_Window)),]#used for significant interactions plot


############ Merge and make plot ###############
MACPET_MACS_Sign_Int_MANGO_data_A5=rbind(MACPET_Sign_Int_MANGO_data_A5,MACS_0.05_Sign_Int_MANGO_data_A5,MACS_0.01_Sign_Int_MANGO_data_A5,MACS_N_Sign_Int_MANGO_data_A5)
save("MACPET_MACS_Sign_Int_MANGO_data_A5",file=file.path(DIRsaveA5,"MACPET_MACS_Sign_Int_MANGO_data_A5"))

# ------------
# significant interactions plot (in row 1):
# ------------
SigInteractions_MANGO_plot_A5=ggplot(MACPET_MACS_Sign_Int_MANGO_data_A5,aes(Merge_Window,Total_significant_int,color=factor(Which),linetype=factor(Which)))+
    geom_line(size=0.3)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.x=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.005, 0),breaks=seq(0,1000,100))+
    ggplot2::scale_y_continuous(expand = c(0.01,0),breaks=seq(min(MACPET_MACS_Sign_Int_MANGO_data_A5$Total_significant_int),
                                                              max(MACPET_MACS_Sign_Int_MANGO_data_A5$Total_significant_int),30))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.9,.35),legend.title=element_blank(),
          legend.text=element_text(size=6),axis.title.x=element_blank(),axis.title.y=element_blank())+guides(col = guide_legend(nrow=2))+
    scale_colour_manual(values=c("red","blue","green","brown"))#coloring factors

SigInteractions_MANGO_plot_A5
save("SigInteractions_MANGO_plot_A5",file=file.path(DIRsaveA5,"SigInteractions_MANGO_plot_A5"))


# ------------
# plots for percentage of peaks inn and peaks used in the interactions(in row 2)
# ------------
PeaksPercentegeInvolved_MANGO_plot_A5=ggplot(MACPET_MACS_Sign_Int_MANGO_data_A5,aes(Merge_Window,100*Peaks_involved_1/Peaks_in,color=factor(Which),linetype=factor(Which)))+
    geom_line(size=0.3)+xlab("MANGO window in bp")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"),
         axis.title.y=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.005, 0),breaks=seq(0,1000,100))+
    ggplot2::scale_y_continuous(expand = c(0.01,0),breaks=seq(min(round(100*MACPET_MACS_Sign_Int_MANGO_data_A5$Peaks_involved_1/MACPET_MACS_Sign_Int_MANGO_data_A5$Peaks_in)),
                                                              max(round(100*MACPET_MACS_Sign_Int_MANGO_data_A5$Peaks_involved_1/MACPET_MACS_Sign_Int_MANGO_data_A5$Peaks_in)),0.2))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.25),legend.title=element_blank(),
          legend.text=element_text(size=6))+guides(col = guide_legend(nrow=2))+
    scale_colour_manual(values=c("red","blue","green","brown"))#coloring factors

PeaksPercentegeInvolved_MANGO_plot_A5
save("PeaksPercentegeInvolved_MANGO_plot_A5",file=file.path(DIRsaveA5,"PeaksPercentegeInvolved_MANGO_plot_A5"))


# ------------
# Plot for interactions overlap for 500 window:
# ------------
# find interactions overlap:
# make Granges MACPET:
Anchor1_MACPET_MANGO_A5=GenomicRanges::GRanges(seqnames=MACPETInteractions500_MANGO_A5$chr_from,
                                               ranges=IRanges::IRanges(start=MACPETInteractions500_MANGO_A5$start_from,
                                                                       end=MACPETInteractions500_MANGO_A5$end_from))
Anchor2_MACPET_MANGO_A5=GenomicRanges::GRanges(seqnames=MACPETInteractions500_MANGO_A5$chr_to,
                                               ranges=IRanges::IRanges(start=MACPETInteractions500_MANGO_A5$start_to,
                                                                       end=MACPETInteractions500_MANGO_A5$end_to))
# make interactions
GInteractions_MACPET_MANGO_A5=GInteractions(Anchor1_MACPET_MANGO_A5,Anchor2_MACPET_MANGO_A5)
# make Granges MACS:
Anchor1_MACS_MANGO_A5=GenomicRanges::GRanges(seqnames=MACSInteractions500_MANGO_A5$chr_from,
                                             ranges=IRanges::IRanges(start=MACSInteractions500_MANGO_A5$start_from,
                                                                     end=MACSInteractions500_MANGO_A5$end_from))
Anchor2_MACS_MANGO_A5=GenomicRanges::GRanges(seqnames=MACSInteractions500_MANGO_A5$chr_to,
                                             ranges=IRanges::IRanges(start=MACSInteractions500_MANGO_A5$start_to,
                                                                     end=MACSInteractions500_MANGO_A5$end_to))

GInteractions_MACS_MANGO_A5=GInteractions(Anchor1_MACS_MANGO_A5,Anchor2_MACS_MANGO_A5)


# make calculations for creating the venn data:
library(ggforce)
TotInteractions_MACPET_MANGO_A5=length(GInteractions_MACPET_MANGO_A5)
TotInteractions_MACS_MANGO_A5=length(GInteractions_MACS_MANGO_A5)
TotInteractions_MANGO_A5=TotInteractions_MACPET_MANGO_A5+TotInteractions_MACS_MANGO_A5
TotCommonInteractions_MACStoMACPET_MANGO_A5=length(InteractionSet::findOverlaps(GInteractions_MACPET_MANGO_A5,GInteractions_MACS_MANGO_A5))
# change to 100 scale
TotInteractions_MACPET_MANGO_A5_100=TotInteractions_MACPET_MANGO_A5*100/TotInteractions_MANGO_A5
TotInteractions_MACS_MANGO_A5_100=TotInteractions_MACS_MANGO_A5*100/TotInteractions_MANGO_A5
TotCommonInteractions_MACStoMACPET_MANGO_A5_100=TotCommonInteractions_MACStoMACPET_MANGO_A5*100/TotInteractions_MANGO_A5

VennInteractions_MANGO_A5_data=data.frame(x0=c(0+TotInteractions_MACPET_MANGO_A5_100/2+TotCommonInteractions_MACStoMACPET_MANGO_A5_100/2,#macpet
                                               100-TotInteractions_MACS_MANGO_A5_100/2-TotCommonInteractions_MACStoMACPET_MANGO_A5_100/2,NA),#macs
                                          y0=c(0.5,0.5,NA),
                                          r=c(TotInteractions_MACPET_MANGO_A5_100/2,TotInteractions_MACS_MANGO_A5_100/2,NA),
                                          Totals=c(TotInteractions_MACPET_MANGO_A5,TotInteractions_MACS_MANGO_A5,TotCommonInteractions_MACStoMACPET_MANGO_A5),
                                          color=c("red","blue",NA),which=c("MACPET","MACS",NA))

save("VennInteractions_MANGO_A5_data",file=file.path(DIRsaveA5,"VennInteractions_MANGO_A5_data"))


VennInteractions_Overlap_MANGO_A5_plot=ggplot()+geom_circle(data=VennInteractions_MANGO_A5_data[2,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="blue",alpha=0.6,size=0.2)+
    geom_circle(data=VennInteractions_MANGO_A5_data[1,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="red",alpha=0.6,size=0.2)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),axis.ticks=element_blank(),
          axis.text=element_blank(),axis.title=element_blank())+
    # add numbers for overlaps:
    annotate("text", x=62, y=0.5, label= as.character(VennInteractions_MANGO_A5_data$Totals[3]),size=2.5)+#common
    annotate("text", x=100-VennInteractions_MANGO_A5_data$r[2], y=0.5, label= as.character(VennInteractions_MANGO_A5_data$Totals[2]-VennInteractions_MANGO_A5_data$Totals[3]),size=2.5)+#macs
    annotate("text", x=0+VennInteractions_MANGO_A5_data$r[1], y=0.5, label= as.character(VennInteractions_MANGO_A5_data$Totals[1]-VennInteractions_MANGO_A5_data$Totals[3]),size=2.5)+#MACPET
    # add groups:
    annotate("text", x=20, y=40, label= "MACPET",size=3)+#MACPET
    annotate("text", x=75, y=40, label= "MACS",size=3)#MACPET


VennInteractions_Overlap_MANGO_A5_plot
save("VennInteractions_Overlap_MANGO_A5_plot",file=file.path(DIRsaveA5,"VennInteractions_Overlap_MANGO_A5_plot"))

# ------------
# Plot for interactions length for 500 window:
# ------------
scientific_10 <- function(x) {
    y=parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
    if(x[1]==0) y[1]=0

    return(y)
}#for scientific values
# Find interaction Distances:
Interactions_Dist_MACPET_MANGO_A5=pairdist(GInteractions_MACPET_MANGO_A5)
Interactions_Dist_MACS_MANGO_A5=pairdist(GInteractions_MACS_MANGO_A5)
Interaction_Dist_Merged_MANGO_A5=data.frame(Dist=c(Interactions_Dist_MACPET_MANGO_A5,Interactions_Dist_MACS_MANGO_A5),
                                            Which=c(rep("MACPET",length(Interactions_Dist_MACPET_MANGO_A5)),
                                                    rep("MACS",length(Interactions_Dist_MACS_MANGO_A5))))
save("Interaction_Dist_Merged_MANGO_A5",file=file.path(DIRsaveA5,"Interaction_Dist_Merged_MANGO_A5"))


Interaction_Dist_MANGO_Plot_A5=ggplot(Interaction_Dist_Merged_MANGO_A5,aes(Dist,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(label=scientific_10,expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"),axis.title.y=element_blank())+
    xlab("Interactions Distance")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

Interaction_Dist_MANGO_Plot_A5
save("Interaction_Dist_MANGO_Plot_A5",file=file.path(DIRsaveA5,"Interaction_Dist_MANGO_Plot_A5"))

#-------------------------------------------------------
#----------------A6: ENCSR000BZY
#-------------------------------------------------------
# The saving directory
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZY"
Repetition="REP1"
DIRsaveA6=file.path(SaveDir,"ENCSR000BZY/MANGO_interactions")
if(!dir.exists(DIRsaveA6)) dir.create(DIRsaveA6)
DIRdata=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/InteractionsResults_MANGO_MACPETInitials/",sep="")# the data directory
############ Load data MACPET ###############
MACPET_files=list.files(file.path(DIRdata,"MACPET"),pattern="Interactions_")
MACPET_Sign_Int_MANGO_data_A6=data.frame()#for significant interactions
MACPET_HeatMap_mat_MANGO_A6=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACPET_HeatMap_mat_MANGO_A6)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACPET_files){#CHANGES
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACPET_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACPET",folder,"MACPET_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACPET_slopPeak_folder=MACPET_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACPET_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACPET_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACPET_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACPET_slopPeak_folder$chr,
                                                  ranges=IRanges::IRanges(start=MACPET_slopPeak_folder$start,
                                                                          end=MACPET_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACPET",folder,"MACPET.interactions.fdr.mango"))){
        MACPET_Interactions_folder=utils::read.table(file.path(DIRdata,"MACPET",folder,"MACPET.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACPET_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACPET_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACPET_Interactions_folder_update=data.frame(chr=c(MACPET_Interactions_folder$chr_from,
                                                           MACPET_Interactions_folder$chr_to),
                                                     start=c(MACPET_Interactions_folder$start_from,
                                                             MACPET_Interactions_folder$start_to),
                                                     end=c(MACPET_Interactions_folder$end_from,
                                                           MACPET_Interactions_folder$end_to))
        # make granges
        MACPET_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACPET_Interactions_folder_update$chr,
                                                                 ranges=IRanges::IRanges(start=MACPET_Interactions_folder_update$start,
                                                                                         end=MACPET_Interactions_folder_update$end))

        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACPET_slopPeak_folder,MACPET_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACPET_HeatMap_mat_MANGO_A6)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACPET_HeatMap_mat_MANGO_A6[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACPET",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#else no interacions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACPET",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }

    MACPET_Sign_Int_MANGO_data_A6=rbind(MACPET_Sign_Int_MANGO_data_A6,Interactions_folder)
}
# sort by merging window
MACPET_Sign_Int_MANGO_data_A6=MACPET_Sign_Int_MANGO_data_A6[with(MACPET_Sign_Int_MANGO_data_A6,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
MACPET_HeatMap_mat_MANGO_A6=MACPET_HeatMap_mat_MANGO_A6[c(1:(Top_max_used)),]#remove rows
MACPET_HeatMap_mat_MANGO_A6[which(is.na(MACPET_HeatMap_mat_MANGO_A6),arr.ind=TRUE)]=0
save("MACPET_HeatMap_mat_MANGO_A6",file=file.path(DIRsaveA6,"MACPET_HeatMap_mat_MANGO_A6"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
MACPETInteractions500_MANGO_A6=utils::read.table(file.path(DIRdata,"MACPET","Interactions_500","MACPET.interactions.fdr.mango"),stringsAsFactors=FALSE)
colnames(MACPETInteractions500_MANGO_A6)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
save("MACPETInteractions500_MANGO_A6",file=file.path(DIRsaveA6,"MACPETInteractions500_MANGO_A6"))

############ Load data MACS ###############
MACS_files=list.files(file.path(DIRdata,"MACS"),pattern="Interactions_")
MACS_0.05_Sign_Int_MANGO_data_A6=data.frame()
MACS_HeatMap_mat_MANGO_A6=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A6)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A6)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A6[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }

        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.05",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.05",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }

    MACS_0.05_Sign_Int_MANGO_data_A6=rbind(MACS_0.05_Sign_Int_MANGO_data_A6,Interactions_folder)
}
# sort by merging window
MACS_0.05_Sign_Int_MANGO_data_A6=MACS_0.05_Sign_Int_MANGO_data_A6[with(MACS_0.05_Sign_Int_MANGO_data_A6,order(Merge_Window)),]#used for significant interactions plot
# reduce matrix:
MACS_HeatMap_mat_MANGO_A6=MACS_HeatMap_mat_MANGO_A6[c(1:(Top_max_used)),]#remove rows
MACS_HeatMap_mat_MANGO_A6[which(is.na(MACS_HeatMap_mat_MANGO_A6),arr.ind=TRUE)]=0
save("MACS_HeatMap_mat_MANGO_A6",file=file.path(DIRsaveA6,"MACS_HeatMap_mat_MANGO_A6"))

# load interactions 500 as suggested window from MANGO to investigate the rest(loads the significant only):
MACSInteractions500_MANGO_A6=utils::read.table(file.path(DIRdata,"MACS","Interactions_500","MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
colnames(MACSInteractions500_MANGO_A6)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
save("MACSInteractions500_MANGO_A6",file=file.path(DIRsaveA6,"MACSInteractions500_MANGO_A6"))

############ Load data MACS ###############
MACS_0.01_files=list.files(file.path(DIRdata,"MACS_0.01"),pattern="Interactions_")
MACS_0.01_Sign_Int_MANGO_data_A6=data.frame()
MACS_HeatMap_mat_MANGO_A6=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A6)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_0.01_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS_0.01",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS_0.01",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS_0.01",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A6)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A6[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.01",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS 0.01",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }
    
    MACS_0.01_Sign_Int_MANGO_data_A6=rbind(MACS_0.01_Sign_Int_MANGO_data_A6,Interactions_folder)
}
# sort by merging window
MACS_0.01_Sign_Int_MANGO_data_A6=MACS_0.01_Sign_Int_MANGO_data_A6[with(MACS_0.01_Sign_Int_MANGO_data_A6,order(Merge_Window)),]#used for significant interactions plot

############ Load data MACS N ###############
MACS_N_files=list.files(file.path(DIRdata,"MACS_N"),pattern="Interactions_")
MACS_N_Sign_Int_MANGO_data_A6=data.frame()
MACS_HeatMap_mat_MANGO_A6=matrix(NA,nrow=100,ncol=11)#for the heat maps
colnames(MACS_HeatMap_mat_MANGO_A6)=seq(0,1000,100)
Top_max_used=1#to remove columns from matrix
for(folder in MACS_N_files){
    # take window:
    Merge_Window=as.numeric(unlist(strsplit(folder,"_"))[2])
    # -------------------
    # load file for the peaks used after merging them:
    # -------------------
    #--------------peaks:
    MACS_slopPeak_folder=utils::read.table(file.path(DIRdata,"MACS_N",folder,"MACS_peaks.slopPeak"),stringsAsFactors=FALSE)
    MACS_slopPeak_folder=MACS_slopPeak_folder[c("V1","V2","V3","V5")]
    colnames(MACS_slopPeak_folder)=c("chr","start","end","speak")
    Peaks_in=nrow(MACS_slopPeak_folder)#peaks in the algorithm
    # make granges
    MACS_slopPeak_folder=GenomicRanges::GRanges(seqnames=MACS_slopPeak_folder$chr,
                                                ranges=IRanges::IRanges(start=MACS_slopPeak_folder$start,
                                                                        end=MACS_slopPeak_folder$end))
    #--------------interactions:
    if(file.exists(file.path(DIRdata,"MACS_N",folder,"MACS.interactions.fdr.mango"))){
        MACS_Interactions_folder=utils::read.table(file.path(DIRdata,"MACS_N",folder,"MACS.interactions.fdr.mango"),stringsAsFactors=FALSE)
        colnames(MACS_Interactions_folder)=c("chr_from","start_from","end_from","chr_to","start_to","end_to","PETs","FDR")
        Total_significant_int=nrow(MACS_Interactions_folder)#significant interactions
        # make one way interactions to find overlaps with the peaks:
        MACS_Interactions_folder_update=data.frame(chr=c(MACS_Interactions_folder$chr_from,
                                                         MACS_Interactions_folder$chr_to),
                                                   start=c(MACS_Interactions_folder$start_from,
                                                           MACS_Interactions_folder$start_to),
                                                   end=c(MACS_Interactions_folder$end_from,
                                                         MACS_Interactions_folder$end_to))
        # make granges
        MACS_Interactions_folder_update=GenomicRanges::GRanges(seqnames=MACS_Interactions_folder_update$chr,
                                                               ranges=IRanges::IRanges(start=MACS_Interactions_folder_update$start,
                                                                                       end=MACS_Interactions_folder_update$end))
        # find overlaps from peaks to interactions:
        Peaks_Used=InteractionSet::countOverlaps(MACS_slopPeak_folder,MACS_Interactions_folder_update)#peaks assosiated with interactions
        # update the matrix:
        Mat_col_index=which(colnames(MACS_HeatMap_mat_MANGO_A6)==Merge_Window)
        max_used=max(Peaks_Used)
        Top_max_used=max(Top_max_used,max_used)
        # fill rows
        for(i in 1:max_used){
            MACS_HeatMap_mat_MANGO_A6[i,Mat_col_index]=length(which(Peaks_Used>=i))
        }
        
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS N",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=Total_significant_int,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=length(which(Peaks_Used>=1)))
    }else{#no interactions
        # create the final data frame:
        Interactions_folder=data.frame(Which="MACS N",
                                       Merge_Window=Merge_Window,
                                       Total_significant_int=0,
                                       Peaks_in=Peaks_in,
                                       Peaks_involved_1=0)
    }
    
    MACS_N_Sign_Int_MANGO_data_A6=rbind(MACS_N_Sign_Int_MANGO_data_A6,Interactions_folder)
}
# sort by merging window
MACS_N_Sign_Int_MANGO_data_A6=MACS_N_Sign_Int_MANGO_data_A6[with(MACS_N_Sign_Int_MANGO_data_A6,order(Merge_Window)),]#used for significant interactions plot


############ Merge and make plot ###############
MACPET_MACS_Sign_Int_MANGO_data_A6=rbind(MACPET_Sign_Int_MANGO_data_A6,MACS_0.05_Sign_Int_MANGO_data_A6,MACS_0.01_Sign_Int_MANGO_data_A6,MACS_N_Sign_Int_MANGO_data_A6)
save("MACPET_MACS_Sign_Int_MANGO_data_A6",file=file.path(DIRsaveA6,"MACPET_MACS_Sign_Int_MANGO_data_A6"))

# ------------
# significant interactions plot (in row 1):
# ------------
SigInteractions_MANGO_plot_A6=ggplot(MACPET_MACS_Sign_Int_MANGO_data_A6,aes(Merge_Window,Total_significant_int,color=factor(Which),linetype=factor(Which)))+
    geom_line(size=0.3)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.x=element_blank(),axis.title.y=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.005, 0),breaks=seq(0,1000,100))+
    ggplot2::scale_y_continuous(expand = c(0.01,0),breaks=seq(min(MACPET_MACS_Sign_Int_MANGO_data_A6$Total_significant_int),
                                                              max(MACPET_MACS_Sign_Int_MANGO_data_A6$Total_significant_int),50))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.35),legend.title=element_blank(),
          legend.text=element_text(size=6))+guides(col = guide_legend(nrow=2))+
    scale_colour_manual(values=c("red","blue","green","brown"))#coloring factors

SigInteractions_MANGO_plot_A6
save("SigInteractions_MANGO_plot_A6",file=file.path(DIRsaveA6,"SigInteractions_MANGO_plot_A6"))


# ------------
# plots for percentage of peaks inn and peaks used in the interactions(in row 2)
# ------------
PeaksPercentegeInvolved_MANGO_plot_A6=ggplot(MACPET_MACS_Sign_Int_MANGO_data_A6,aes(Merge_Window,100*Peaks_involved_1/Peaks_in,color=factor(Which),linetype=factor(Which)))+
    geom_line(size=0.3)+xlab(" ")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"),axis.title.y=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.005, 0),breaks=seq(0,1000,100))+
    ggplot2::scale_y_continuous(expand = c(0.01,0),breaks=seq(min(round(100*MACPET_MACS_Sign_Int_MANGO_data_A6$Peaks_involved_1/MACPET_MACS_Sign_Int_MANGO_data_A6$Peaks_in)),
                                                              max(round(100*MACPET_MACS_Sign_Int_MANGO_data_A6$Peaks_involved_1/MACPET_MACS_Sign_Int_MANGO_data_A6$Peaks_in)),1))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.35),legend.title=element_blank(),
          legend.text=element_text(size=6))+guides(col = guide_legend(nrow=2))+
    scale_colour_manual(values=c("red","blue","green","brown"))#coloring factors

PeaksPercentegeInvolved_MANGO_plot_A6
save("PeaksPercentegeInvolved_MANGO_plot_A6",file=file.path(DIRsaveA6,"PeaksPercentegeInvolved_MANGO_plot_A6"))

# ------------
# Plot for interactions overlap for 500 window:
# ------------
# find interactions overlap:
# make Granges MACPET:
Anchor1_MACPET_MANGO_A6=GenomicRanges::GRanges(seqnames=MACPETInteractions500_MANGO_A6$chr_from,
                                               ranges=IRanges::IRanges(start=MACPETInteractions500_MANGO_A6$start_from,
                                                                       end=MACPETInteractions500_MANGO_A6$end_from))
Anchor2_MACPET_MANGO_A6=GenomicRanges::GRanges(seqnames=MACPETInteractions500_MANGO_A6$chr_to,
                                               ranges=IRanges::IRanges(start=MACPETInteractions500_MANGO_A6$start_to,
                                                                       end=MACPETInteractions500_MANGO_A6$end_to))
# make interactions
GInteractions_MACPET_MANGO_A6=GInteractions(Anchor1_MACPET_MANGO_A6,Anchor2_MACPET_MANGO_A6)
# make Granges MACS:
Anchor1_MACS_MANGO_A6=GenomicRanges::GRanges(seqnames=MACSInteractions500_MANGO_A6$chr_from,
                                             ranges=IRanges::IRanges(start=MACSInteractions500_MANGO_A6$start_from,
                                                                     end=MACSInteractions500_MANGO_A6$end_from))
Anchor2_MACS_MANGO_A6=GenomicRanges::GRanges(seqnames=MACSInteractions500_MANGO_A6$chr_to,
                                             ranges=IRanges::IRanges(start=MACSInteractions500_MANGO_A6$start_to,
                                                                     end=MACSInteractions500_MANGO_A6$end_to))

GInteractions_MACS_MANGO_A6=GInteractions(Anchor1_MACS_MANGO_A6,Anchor2_MACS_MANGO_A6)


# make calculations for creating the venn data:
library(ggforce)
TotInteractions_MACPET_MANGO_A6=length(GInteractions_MACPET_MANGO_A6)
TotInteractions_MACS_MANGO_A6=length(GInteractions_MACS_MANGO_A6)
TotInteractions_MANGO_A6=TotInteractions_MACPET_MANGO_A6+TotInteractions_MACS_MANGO_A6
TotCommonInteractions_MACStoMACPET_MANGO_A6=length(InteractionSet::findOverlaps(GInteractions_MACPET_MANGO_A6,GInteractions_MACS_MANGO_A6))
# change to 100 scale
TotInteractions_MACPET_MANGO_A6_100=TotInteractions_MACPET_MANGO_A6*100/TotInteractions_MANGO_A6
TotInteractions_MACS_MANGO_A6_100=TotInteractions_MACS_MANGO_A6*100/TotInteractions_MANGO_A6
TotCommonInteractions_MACStoMACPET_MANGO_A6_100=TotCommonInteractions_MACStoMACPET_MANGO_A6*100/TotInteractions_MANGO_A6

VennInteractions_MANGO_A6_data=data.frame(x0=c(0+TotInteractions_MACPET_MANGO_A6_100/2+TotCommonInteractions_MACStoMACPET_MANGO_A6_100/2,#macpet
                                               100-TotInteractions_MACS_MANGO_A6_100/2-TotCommonInteractions_MACStoMACPET_MANGO_A6_100/2,NA),#macs
                                          y0=c(0.5,0.5,NA),
                                          r=c(TotInteractions_MACPET_MANGO_A6_100/2,TotInteractions_MACS_MANGO_A6_100/2,NA),
                                          Totals=c(TotInteractions_MACPET_MANGO_A6,TotInteractions_MACS_MANGO_A6,TotCommonInteractions_MACStoMACPET_MANGO_A6),
                                          color=c("red","blue",NA),which=c("MACPET","MACS",NA))

save("VennInteractions_MANGO_A6_data",file=file.path(DIRsaveA6,"VennInteractions_MANGO_A6_data"))


VennInteractions_Overlap_MANGO_A6_plot=ggplot()+geom_circle(data=VennInteractions_MANGO_A6_data[2,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="blue",alpha=0.6,size=0.2)+
    geom_circle(data=VennInteractions_MANGO_A6_data[1,],mapping=aes(x0=x0,y0=y0,r=r),color="black",fill="red",alpha=0.6,size=0.2)+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),axis.ticks=element_blank(),
          axis.text=element_blank(),axis.title=element_blank())+
    # add numbers for overlaps:
    annotate("text", x=77, y=0.5, label= as.character(VennInteractions_MANGO_A6_data$Totals[3]),size=2.5)+#common
    annotate("text", x=100-VennInteractions_MANGO_A6_data$r[2], y=0.5, label= as.character(VennInteractions_MANGO_A6_data$Totals[2]-VennInteractions_MANGO_A6_data$Totals[3]),size=2.5)+#macs
    annotate("text", x=0+VennInteractions_MANGO_A6_data$r[1], y=0.5, label= as.character(VennInteractions_MANGO_A6_data$Totals[1]-VennInteractions_MANGO_A6_data$Totals[3]),size=2.5)+#MACPET
    # add groups:
    annotate("text", x=15, y=40, label= "MACPET",size=3)+#MACPET
    annotate("text", x=80, y=40, label= "MACS",size=3)#MACPET


VennInteractions_Overlap_MANGO_A6_plot
save("VennInteractions_Overlap_MANGO_A6_plot",file=file.path(DIRsaveA6,"VennInteractions_Overlap_MANGO_A6_plot"))

# ------------
# Plot for interactions length for 500 window:
# ------------
scientific_10 <- function(x) {
    y=parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
    if(x[1]==0) y[1]=0

    return(y)
}#for scientific values
# Find interaction Distances:
Interactions_Dist_MACPET_MANGO_A6=pairdist(GInteractions_MACPET_MANGO_A6)
Interactions_Dist_MACS_MANGO_A6=pairdist(GInteractions_MACS_MANGO_A6)
Interaction_Dist_Merged_MANGO_A6=data.frame(Dist=c(Interactions_Dist_MACPET_MANGO_A6,Interactions_Dist_MACS_MANGO_A6),
                                            Which=c(rep("MACPET",length(Interactions_Dist_MACPET_MANGO_A6)),
                                                    rep("MACS",length(Interactions_Dist_MACS_MANGO_A6))))
save("Interaction_Dist_Merged_MANGO_A6",file=file.path(DIRsaveA6,"Interaction_Dist_Merged_MANGO_A6"))


Interaction_Dist_MANGO_Plot_A6=ggplot(Interaction_Dist_Merged_MANGO_A6,aes(Dist,color=factor(Which),linetype=factor(Which),fill=factor(Which)))+geom_density(alpha=0.6)+
    # scale of y and removing space beneath:
    ggplot2::scale_y_continuous(label=scientific_10,expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0.001, 0))+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=6,family="Times"),
          axis.title=element_text(size=10,family="Times"),axis.title.y=element_blank())+
    xlab("")+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.75),legend.title=element_blank(),
          legend.text=element_text(size=7))+
    scale_colour_manual(values=c("red","blue"))#coloring factors

Interaction_Dist_MANGO_Plot_A6
save("Interaction_Dist_MANGO_Plot_A6",file=file.path(DIRsaveA6,"Interaction_Dist_MANGO_Plot_A6"))


#################################################################################################################
################################ Figures for MO, SE, FDR for ALL the parameters: ######################################
#################################################################################################################
# The saving directory
#----------------
#----------------A1: ENCSR000BZZ
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
############ Load data MACPET ###############
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000BZZ"
Repetition="REP1"
DIRsaveA1=file.path(SaveDir,"ENCSR000BZZ/MultipleMotvComp")
if(!dir.exists(DIRsaveA1)) dir.create(DIRsaveA1)
#save in the external drive:
DIR=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results",sep="")# the data directory
setwd(DIR)
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACPET.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACPET")
    data.MACPET.run_i=data.motif.overlap.expected
    data.MACPET.run_i$Run=i
    data.MACPET.Allruns=rbind(data.MACPET.Allruns,data.MACPET.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACPET.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACPET.run_1=subset(data.MACPET.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACPET.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACPET.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACPET.run_2=subset(data.MACPET.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACPET.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACPET.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACPET.run_3=subset(data.MACPET.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACPET.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACPET.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACPET.run_4=subset(data.MACPET.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACPET.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACPET.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACPET.run_5=subset(data.MACPET.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACPET.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACPET.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACPET.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                      Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                      Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                      Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                      Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                      Spatial.error.run_1=Spatial.error.run_1,
                                                      Spatial.error.run_2=Spatial.error.run_2,
                                                      Spatial.error.run_3=Spatial.error.run_3,
                                                      Spatial.error.run_4=Spatial.error.run_4,
                                                      Spatial.error.run_5=Spatial.error.run_5,
                                                      Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                      Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                      Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                  tot.motifs.run_2/tot.bs,
                                                                                  tot.motifs.run_3/tot.bs,
                                                                                  tot.motifs.run_4/tot.bs,
                                                                                  tot.motifs.run_5/tot.bs)),
                                                      Spatial.error.min=min(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                      Spatial.error.max=max(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                      Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                                Spatial.error.run_2,
                                                                                Spatial.error.run_3,
                                                                                Spatial.error.run_4,
                                                                                Spatial.error.run_5)),
                                                      tot.bs=x,Which="MACPET")
    
    # merge:
    data.MACPET.occurrence.spatial.error=rbind(data.MACPET.occurrence.spatial.error,
                                               data.MACPET.occurrence.spatial.error.x)
    
}
setwd("../../../..")
getwd()
############ Load data MACS for default###############
setwd("MACS_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS.run_i=data.motif.overlap.expected
    data.MACS.run_i$Run=i
    data.MACS.Allruns=rbind(data.MACS.Allruns,data.MACS.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                    Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                    Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                    Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                    Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                    Spatial.error.run_1=Spatial.error.run_1,
                                                    Spatial.error.run_2=Spatial.error.run_2,
                                                    Spatial.error.run_3=Spatial.error.run_3,
                                                    Spatial.error.run_4=Spatial.error.run_4,
                                                    Spatial.error.run_5=Spatial.error.run_5,
                                                    Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                              tot.motifs.run_2/tot.bs,
                                                                              tot.motifs.run_3/tot.bs,
                                                                              tot.motifs.run_4/tot.bs,
                                                                              tot.motifs.run_5/tot.bs)),
                                                    Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                              tot.motifs.run_2/tot.bs,
                                                                              tot.motifs.run_3/tot.bs,
                                                                              tot.motifs.run_4/tot.bs,
                                                                              tot.motifs.run_5/tot.bs)),
                                                    Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                    Spatial.error.min=min(c(Spatial.error.run_1,
                                                                            Spatial.error.run_2,
                                                                            Spatial.error.run_3,
                                                                            Spatial.error.run_4,
                                                                            Spatial.error.run_5)),
                                                    Spatial.error.max=max(c(Spatial.error.run_1,
                                                                            Spatial.error.run_2,
                                                                            Spatial.error.run_3,
                                                                            Spatial.error.run_4,
                                                                            Spatial.error.run_5)),
                                                    Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                    tot.bs=x,Which="MACS")
    
    # merge:
    data.MACS.occurrence.spatial.error=rbind(data.MACS.occurrence.spatial.error,
                                             data.MACS.occurrence.spatial.error.x)
    
}
setwd("../../..")
getwd()
############ Load data MACS2 for default###############
setwd("MACS2_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS2.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS2.run_i=data.motif.overlap.expected
    data.MACS2.run_i$Run=i
    data.MACS2.Allruns=rbind(data.MACS2.Allruns,data.MACS2.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS2.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS2.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS2.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS2.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS2.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS2.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS2.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                    Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                    Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                    Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                    Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                    Spatial.error.run_1=Spatial.error.run_1,
                                                    Spatial.error.run_2=Spatial.error.run_2,
                                                    Spatial.error.run_3=Spatial.error.run_3,
                                                    Spatial.error.run_4=Spatial.error.run_4,
                                                    Spatial.error.run_5=Spatial.error.run_5,
                                                    Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                              tot.motifs.run_2/tot.bs,
                                                                              tot.motifs.run_3/tot.bs,
                                                                              tot.motifs.run_4/tot.bs,
                                                                              tot.motifs.run_5/tot.bs)),
                                                    Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                              tot.motifs.run_2/tot.bs,
                                                                              tot.motifs.run_3/tot.bs,
                                                                              tot.motifs.run_4/tot.bs,
                                                                              tot.motifs.run_5/tot.bs)),
                                                    Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                    Spatial.error.min=min(c(Spatial.error.run_1,
                                                                            Spatial.error.run_2,
                                                                            Spatial.error.run_3,
                                                                            Spatial.error.run_4,
                                                                            Spatial.error.run_5)),
                                                    Spatial.error.max=max(c(Spatial.error.run_1,
                                                                            Spatial.error.run_2,
                                                                            Spatial.error.run_3,
                                                                            Spatial.error.run_4,
                                                                            Spatial.error.run_5)),
                                                    Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                    tot.bs=x,Which="MACS(2)")
    
    # merge:
    data.MACS2.occurrence.spatial.error=rbind(data.MACS2.occurrence.spatial.error,
                                             data.MACS2.occurrence.spatial.error.x)
    
}
setwd("../../..")
getwd()
############ Load data MACS3 for default###############
setwd("MACS3_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS3.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS3.run_i=data.motif.overlap.expected
    data.MACS3.run_i$Run=i
    data.MACS3.Allruns=rbind(data.MACS3.Allruns,data.MACS3.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS3.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS3.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS3.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS3.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS3.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS3.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS3.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                     Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                     Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                     Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                     Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                     Spatial.error.run_1=Spatial.error.run_1,
                                                     Spatial.error.run_2=Spatial.error.run_2,
                                                     Spatial.error.run_3=Spatial.error.run_3,
                                                     Spatial.error.run_4=Spatial.error.run_4,
                                                     Spatial.error.run_5=Spatial.error.run_5,
                                                     Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                 tot.motifs.run_2/tot.bs,
                                                                                 tot.motifs.run_3/tot.bs,
                                                                                 tot.motifs.run_4/tot.bs,
                                                                                 tot.motifs.run_5/tot.bs)),
                                                     Spatial.error.min=min(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.max=max(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                               Spatial.error.run_2,
                                                                               Spatial.error.run_3,
                                                                               Spatial.error.run_4,
                                                                               Spatial.error.run_5)),
                                                     tot.bs=x,Which="MACS(3)")
    
    # merge:
    data.MACS3.occurrence.spatial.error=rbind(data.MACS3.occurrence.spatial.error,
                                              data.MACS3.occurrence.spatial.error.x)
    
}
setwd("../../..")
getwd()
############ Load data MACS4 for default###############
setwd("MACS4_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS4.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS4.run_i=data.motif.overlap.expected
    data.MACS4.run_i$Run=i
    data.MACS4.Allruns=rbind(data.MACS4.Allruns,data.MACS4.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS4.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS4.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS4.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS4.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS4.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS4.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS4.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                     Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                     Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                     Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                     Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                     Spatial.error.run_1=Spatial.error.run_1,
                                                     Spatial.error.run_2=Spatial.error.run_2,
                                                     Spatial.error.run_3=Spatial.error.run_3,
                                                     Spatial.error.run_4=Spatial.error.run_4,
                                                     Spatial.error.run_5=Spatial.error.run_5,
                                                     Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                 tot.motifs.run_2/tot.bs,
                                                                                 tot.motifs.run_3/tot.bs,
                                                                                 tot.motifs.run_4/tot.bs,
                                                                                 tot.motifs.run_5/tot.bs)),
                                                     Spatial.error.min=min(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.max=max(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                               Spatial.error.run_2,
                                                                               Spatial.error.run_3,
                                                                               Spatial.error.run_4,
                                                                               Spatial.error.run_5)),
                                                     tot.bs=x,Which="MACS(4)")
    
    # merge:
    data.MACS4.occurrence.spatial.error=rbind(data.MACS4.occurrence.spatial.error,
                                              data.MACS4.occurrence.spatial.error.x)
    
}
setwd("../../..")
getwd()
############ Create the whole data ###############
Motif.occurence.allA1=rbind(data.MACPET.occurrence.spatial.error,data.MACS.occurrence.spatial.error,
                            data.MACS2.occurrence.spatial.error,
                            data.MACS3.occurrence.spatial.error,
                            data.MACS4.occurrence.spatial.error)
save("Motif.occurence.allA1",file=file.path(DIRsaveA1,"Motif.occurence.allA1"))
############ MO plot ###############

MOggplotA1=ggplot(Motif.occurence.allA1,aes(tot.bs,100*Motif.occurence.mean,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    ylab("Motif Occurance (%)")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.x=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0, 0))+ggplot2::scale_y_continuous(expand = c(0, 0),breaks=seq(from=0,to=100,by=5),limits=c(55,85))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.50),legend.title=element_blank(),
          legend.text=element_text(size=6))+guides(col = guide_legend(ncol=2))+
    scale_colour_manual(values=c("red","blue","green","black","purple"))#coloring factors


# MOggplotA1
save("MOggplotA1",file=file.path(DIRsaveA1,"MOggplotA1"))

############ SE plot ###############
SEggplotA1=ggplot(Motif.occurence.allA1,aes(tot.bs,Spatial.error.mean,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    ylab("Spatial Resolution (bp)")+xlab(" ")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    ggplot2::scale_x_continuous(expand = c(0, 0))+ggplot2::scale_y_continuous(expand = c(0, 0),breaks=seq(from=0,to=100,by=5),limits=c(15,35))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.85,.0),legend.title=element_blank(),
          legend.text=element_text(size=6))+guides(col = guide_legend(ncol=2))+
    scale_colour_manual(values=c("red","blue","green","black","purple"))#coloring factors

# SEggplotA1
save("SEggplotA1",file=file.path(DIRsaveA1,"SEggplotA1"))

############ FDR plot ###############
MACPET_run_1=subset(data.MACPET.Allruns,Run==1)
MACPET_run_1$FDR=MACPET_run_1$FDR*100/max(MACPET_run_1$FDR)
MACS_run_1=subset(data.MACS.Allruns,Run==1)
MACS_run_1$p.values.FDR=MACS_run_1$p.values.FDR*100/max(MACS_run_1$p.values.FDR)
MACS2_run_1=subset(data.MACS2.Allruns,Run==1)
MACS2_run_1$p.values.FDR=MACS2_run_1$p.values.FDR*100/max(MACS2_run_1$p.values.FDR)
MACS3_run_1=subset(data.MACS3.Allruns,Run==1)
MACS3_run_1$p.values.FDR=MACS3_run_1$p.values.FDR*100/max(MACS3_run_1$p.values.FDR)
MACS4_run_1=subset(data.MACS4.Allruns,Run==1)
MACS4_run_1$p.values.FDR=MACS4_run_1$p.values.FDR*100/max(MACS4_run_1$p.values.FDR)

data.FDRA1=data.frame(FDR=c(sort(MACPET_run_1$FDR),sort(MACS_run_1$p.values.FDR),
                            sort(MACS2_run_1$p.values.FDR),sort(MACS3_run_1$p.values.FDR),sort(MACS4_run_1$p.values.FDR)),
                      Which=c(rep("MACPET",nrow(MACPET_run_1)),rep("MACS",nrow(MACS_run_1)),
                              rep("MACS(2)",nrow(MACS2_run_1)),rep("MACS(3)",nrow(MACS3_run_1)),
                              rep("MACS(4)",nrow(MACS4_run_1))),
                      bs=c(1:nrow(MACPET_run_1),1:nrow(MACS_run_1),
                           1:nrow(MACS2_run_1),1:nrow(MACS3_run_1),1:nrow(MACS4_run_1)))

save("data.FDRA1",file=file.path(DIRsaveA1,"data.FDRA1"))

FDRggplotA1=ggplot(data.FDRA1,aes(bs,FDR,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    ylab("FDR (%)")+xlab(" ")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"))+
    ggplot2::scale_x_continuous(expand = c(0, 0))+ggplot2::scale_y_continuous(expand = c(0.005, 0),breaks=seq(from=0,to=100,by=10),limits=c(0,100))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.55,.30),legend.title=element_blank(),
          legend.text=element_text(size=6))+
    scale_colour_manual(values=c("red","blue","green","black","purple"))#coloring factors

# FDRggplotA1
save("FDRggplotA1",file=file.path(DIRsaveA1,"FDRggplotA1"))

#################################################################################################################
################################ Figures for MO, SE, FDR for ALL the parameters: ######################################
#################################################################################################################
# The saving directory
#----------------
#----------------A2: ENCSR000CAD
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
############ Load data MACPET ###############
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAD"
Repetition="REP1"
DIRsaveA1=file.path(SaveDir,"ENCSR000CAD/MultipleMotvComp")
if(!dir.exists(DIRsaveA1)) dir.create(DIRsaveA1)
#save in the external drive:
DIR=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results",sep="")# the data directory
setwd(DIR)
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACPET.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACPET")
    data.MACPET.run_i=data.motif.overlap.expected
    data.MACPET.run_i$Run=i
    data.MACPET.Allruns=rbind(data.MACPET.Allruns,data.MACPET.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACPET.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACPET.run_1=subset(data.MACPET.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACPET.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACPET.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACPET.run_2=subset(data.MACPET.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACPET.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACPET.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACPET.run_3=subset(data.MACPET.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACPET.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACPET.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACPET.run_4=subset(data.MACPET.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACPET.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACPET.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACPET.run_5=subset(data.MACPET.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACPET.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACPET.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACPET.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                      Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                      Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                      Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                      Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                      Spatial.error.run_1=Spatial.error.run_1,
                                                      Spatial.error.run_2=Spatial.error.run_2,
                                                      Spatial.error.run_3=Spatial.error.run_3,
                                                      Spatial.error.run_4=Spatial.error.run_4,
                                                      Spatial.error.run_5=Spatial.error.run_5,
                                                      Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                      Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                      Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                  tot.motifs.run_2/tot.bs,
                                                                                  tot.motifs.run_3/tot.bs,
                                                                                  tot.motifs.run_4/tot.bs,
                                                                                  tot.motifs.run_5/tot.bs)),
                                                      Spatial.error.min=min(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                      Spatial.error.max=max(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                      Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                                Spatial.error.run_2,
                                                                                Spatial.error.run_3,
                                                                                Spatial.error.run_4,
                                                                                Spatial.error.run_5)),
                                                      tot.bs=x,Which="MACPET")
    
    # merge:
    data.MACPET.occurrence.spatial.error=rbind(data.MACPET.occurrence.spatial.error,
                                               data.MACPET.occurrence.spatial.error.x)
    
}
setwd("../../../..")
getwd()
############ Load data MACS for default###############
setwd("MACS_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS.run_i=data.motif.overlap.expected
    data.MACS.run_i$Run=i
    data.MACS.Allruns=rbind(data.MACS.Allruns,data.MACS.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                    Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                    Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                    Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                    Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                    Spatial.error.run_1=Spatial.error.run_1,
                                                    Spatial.error.run_2=Spatial.error.run_2,
                                                    Spatial.error.run_3=Spatial.error.run_3,
                                                    Spatial.error.run_4=Spatial.error.run_4,
                                                    Spatial.error.run_5=Spatial.error.run_5,
                                                    Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                              tot.motifs.run_2/tot.bs,
                                                                              tot.motifs.run_3/tot.bs,
                                                                              tot.motifs.run_4/tot.bs,
                                                                              tot.motifs.run_5/tot.bs)),
                                                    Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                              tot.motifs.run_2/tot.bs,
                                                                              tot.motifs.run_3/tot.bs,
                                                                              tot.motifs.run_4/tot.bs,
                                                                              tot.motifs.run_5/tot.bs)),
                                                    Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                    Spatial.error.min=min(c(Spatial.error.run_1,
                                                                            Spatial.error.run_2,
                                                                            Spatial.error.run_3,
                                                                            Spatial.error.run_4,
                                                                            Spatial.error.run_5)),
                                                    Spatial.error.max=max(c(Spatial.error.run_1,
                                                                            Spatial.error.run_2,
                                                                            Spatial.error.run_3,
                                                                            Spatial.error.run_4,
                                                                            Spatial.error.run_5)),
                                                    Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                    tot.bs=x,Which="MACS")
    
    # merge:
    data.MACS.occurrence.spatial.error=rbind(data.MACS.occurrence.spatial.error,
                                             data.MACS.occurrence.spatial.error.x)
    
}
setwd("../../..")
getwd()
############ Load data MACS2 for default###############
setwd("MACS2_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS2.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS2.run_i=data.motif.overlap.expected
    data.MACS2.run_i$Run=i
    data.MACS2.Allruns=rbind(data.MACS2.Allruns,data.MACS2.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS2.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS2.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS2.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS2.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS2.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS2.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS2.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                     Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                     Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                     Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                     Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                     Spatial.error.run_1=Spatial.error.run_1,
                                                     Spatial.error.run_2=Spatial.error.run_2,
                                                     Spatial.error.run_3=Spatial.error.run_3,
                                                     Spatial.error.run_4=Spatial.error.run_4,
                                                     Spatial.error.run_5=Spatial.error.run_5,
                                                     Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                 tot.motifs.run_2/tot.bs,
                                                                                 tot.motifs.run_3/tot.bs,
                                                                                 tot.motifs.run_4/tot.bs,
                                                                                 tot.motifs.run_5/tot.bs)),
                                                     Spatial.error.min=min(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.max=max(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                               Spatial.error.run_2,
                                                                               Spatial.error.run_3,
                                                                               Spatial.error.run_4,
                                                                               Spatial.error.run_5)),
                                                     tot.bs=x,Which="MACS(2)")
    
    # merge:
    data.MACS2.occurrence.spatial.error=rbind(data.MACS2.occurrence.spatial.error,
                                              data.MACS2.occurrence.spatial.error.x)
    
}
setwd("../../..")
getwd()
############ Load data MACS3 for default###############
setwd("MACS3_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS3.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS3.run_i=data.motif.overlap.expected
    data.MACS3.run_i$Run=i
    data.MACS3.Allruns=rbind(data.MACS3.Allruns,data.MACS3.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS3.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS3.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS3.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS3.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS3.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS3.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS3.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                     Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                     Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                     Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                     Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                     Spatial.error.run_1=Spatial.error.run_1,
                                                     Spatial.error.run_2=Spatial.error.run_2,
                                                     Spatial.error.run_3=Spatial.error.run_3,
                                                     Spatial.error.run_4=Spatial.error.run_4,
                                                     Spatial.error.run_5=Spatial.error.run_5,
                                                     Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                 tot.motifs.run_2/tot.bs,
                                                                                 tot.motifs.run_3/tot.bs,
                                                                                 tot.motifs.run_4/tot.bs,
                                                                                 tot.motifs.run_5/tot.bs)),
                                                     Spatial.error.min=min(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.max=max(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                               Spatial.error.run_2,
                                                                               Spatial.error.run_3,
                                                                               Spatial.error.run_4,
                                                                               Spatial.error.run_5)),
                                                     tot.bs=x,Which="MACS(3)")
    
    # merge:
    data.MACS3.occurrence.spatial.error=rbind(data.MACS3.occurrence.spatial.error,
                                              data.MACS3.occurrence.spatial.error.x)
    
}
setwd("../../..")
getwd()
############ Load data MACS4 for default###############
setwd("MACS4_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS4.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS4.run_i=data.motif.overlap.expected
    data.MACS4.run_i$Run=i
    data.MACS4.Allruns=rbind(data.MACS4.Allruns,data.MACS4.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS4.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS4.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS4.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS4.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS4.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS4.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS4.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                     Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                     Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                     Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                     Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                     Spatial.error.run_1=Spatial.error.run_1,
                                                     Spatial.error.run_2=Spatial.error.run_2,
                                                     Spatial.error.run_3=Spatial.error.run_3,
                                                     Spatial.error.run_4=Spatial.error.run_4,
                                                     Spatial.error.run_5=Spatial.error.run_5,
                                                     Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                 tot.motifs.run_2/tot.bs,
                                                                                 tot.motifs.run_3/tot.bs,
                                                                                 tot.motifs.run_4/tot.bs,
                                                                                 tot.motifs.run_5/tot.bs)),
                                                     Spatial.error.min=min(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.max=max(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                               Spatial.error.run_2,
                                                                               Spatial.error.run_3,
                                                                               Spatial.error.run_4,
                                                                               Spatial.error.run_5)),
                                                     tot.bs=x,Which="MACS(4)")
    
    # merge:
    data.MACS4.occurrence.spatial.error=rbind(data.MACS4.occurrence.spatial.error,
                                              data.MACS4.occurrence.spatial.error.x)
    
}
setwd("../../..")
getwd()
############ Create the whole data ###############
Motif.occurence.allA2=rbind(data.MACPET.occurrence.spatial.error,data.MACS.occurrence.spatial.error,
                            data.MACS2.occurrence.spatial.error,
                            data.MACS3.occurrence.spatial.error,
                            data.MACS4.occurrence.spatial.error)
save("Motif.occurence.allA2",file=file.path(DIRsaveA1,"Motif.occurence.allA2"))
############ MO plot ###############

MOggplotA2=ggplot(Motif.occurence.allA2,aes(tot.bs,100*Motif.occurence.mean,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
   # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.y=element_blank(),
          axis.title.x=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0, 0),breaks=seq(from=0,to=100,by=3),limits=c(80,95))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.0),legend.title=element_blank(),
          legend.text=element_text(size=6))+guides(col = guide_legend(ncol=2))+
    scale_colour_manual(values=c("red","blue","green","black","purple"))#coloring factors


# MOggplotA2
save("MOggplotA2",file=file.path(DIRsaveA1,"MOggplotA2"))

############ SE plot ###############
SEggplotA2=ggplot(Motif.occurence.allA2,aes(tot.bs,Spatial.error.mean,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
# remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.y=element_blank())+xlab("Sorted Binding Sites")+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0, 0),breaks=seq(from=0,to=100,by=3),limits=c(10,30))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.60),legend.title=element_blank(),
          legend.text=element_text(size=6))+guides(col = guide_legend(ncol=2))+
    scale_colour_manual(values=c("red","blue","green","black","purple"))#coloring factors

# SEggplotA2
save("SEggplotA2",file=file.path(DIRsaveA1,"SEggplotA2"))

############ FDR plot ###############
MACPET_run_1=subset(data.MACPET.Allruns,Run==1)
MACPET_run_1$FDR=MACPET_run_1$FDR*100/max(MACPET_run_1$FDR)
MACS_run_1=subset(data.MACS.Allruns,Run==1)
MACS_run_1$p.values.FDR=MACS_run_1$p.values.FDR*100/max(MACS_run_1$p.values.FDR)
MACS2_run_1=subset(data.MACS2.Allruns,Run==1)
MACS2_run_1$p.values.FDR=MACS2_run_1$p.values.FDR*100/max(MACS2_run_1$p.values.FDR)
MACS3_run_1=subset(data.MACS3.Allruns,Run==1)
MACS3_run_1$p.values.FDR=MACS3_run_1$p.values.FDR*100/max(MACS3_run_1$p.values.FDR)
MACS4_run_1=subset(data.MACS4.Allruns,Run==1)
MACS4_run_1$p.values.FDR=MACS4_run_1$p.values.FDR*100/max(MACS4_run_1$p.values.FDR)


data.FDRA2=data.frame(FDR=c(sort(MACPET_run_1$FDR),sort(MACS_run_1$p.values.FDR),
                            sort(MACS2_run_1$p.values.FDR),sort(MACS3_run_1$p.values.FDR),sort(MACS4_run_1$p.values.FDR)),
                      Which=c(rep("MACPET",nrow(MACPET_run_1)),rep("MACS",nrow(MACS_run_1)),
                              rep("MACS(2)",nrow(MACS2_run_1)),rep("MACS(3)",nrow(MACS3_run_1)),
                              rep("MACS(4)",nrow(MACS4_run_1))),
                      bs=c(1:nrow(MACPET_run_1),1:nrow(MACS_run_1),
                           1:nrow(MACS2_run_1),1:nrow(MACS3_run_1),1:nrow(MACS4_run_1)))

save("data.FDRA2",file=file.path(DIRsaveA1,"data.FDRA2"))

FDRggplotA2=ggplot(data.FDRA2,aes(bs,FDR,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    xlab("Sorted Binding Sites")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.y=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0.005, 0),breaks=seq(from=0,to=100,by=10),limits=c(0,100))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.55,.30),legend.title=element_blank(),
          legend.text=element_text(size=6))+
    scale_colour_manual(values=c("red","blue","green","black","purple"))#coloring factors

# FDRggplotA2
save("FDRggplotA2",file=file.path(DIRsaveA1,"FDRggplotA2"))

#################################################################################################################
################################ Figures for MO, SE, FDR for ALL the parameters: ######################################
#################################################################################################################
# The saving directory
#----------------
#----------------A3: ENCSR000CAC
#----------------
SaveDir="A_PATH/Figures"
dir.exists(SaveDir)
############ Load data MACPET ###############
DATA.type="ChIA_PET_Data"
Experiment="Experiment_ENCSR000CAC"
Repetition="REP1"
DIRsaveA1=file.path(SaveDir,"ENCSR000CAC/MultipleMotvComp")
if(!dir.exists(DIRsaveA1)) dir.create(DIRsaveA1)
#save in the external drive:
DIR=paste("YOURPATH/",DATA.type,"/",Experiment,"/",Repetition,"/MACPET_results/S3_results",sep="")# the data directory
setwd(DIR)
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACPET.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACPET")
    data.MACPET.run_i=data.motif.overlap.expected
    data.MACPET.run_i$Run=i
    data.MACPET.Allruns=rbind(data.MACPET.Allruns,data.MACPET.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACPET.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACPET.run_1=subset(data.MACPET.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACPET.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACPET.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACPET.run_2=subset(data.MACPET.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACPET.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACPET.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACPET.run_3=subset(data.MACPET.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACPET.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACPET.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACPET.run_4=subset(data.MACPET.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACPET.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACPET.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACPET.run_5=subset(data.MACPET.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACPET.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACPET.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACPET.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                      Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                      Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                      Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                      Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                      Spatial.error.run_1=Spatial.error.run_1,
                                                      Spatial.error.run_2=Spatial.error.run_2,
                                                      Spatial.error.run_3=Spatial.error.run_3,
                                                      Spatial.error.run_4=Spatial.error.run_4,
                                                      Spatial.error.run_5=Spatial.error.run_5,
                                                      Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                      Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                      Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                  tot.motifs.run_2/tot.bs,
                                                                                  tot.motifs.run_3/tot.bs,
                                                                                  tot.motifs.run_4/tot.bs,
                                                                                  tot.motifs.run_5/tot.bs)),
                                                      Spatial.error.min=min(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                      Spatial.error.max=max(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                      Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                                Spatial.error.run_2,
                                                                                Spatial.error.run_3,
                                                                                Spatial.error.run_4,
                                                                                Spatial.error.run_5)),
                                                      tot.bs=x,Which="MACPET")
    
    # merge:
    data.MACPET.occurrence.spatial.error=rbind(data.MACPET.occurrence.spatial.error,
                                               data.MACPET.occurrence.spatial.error.x)
    
}
setwd("../../../..")
getwd()
############ Load data MACS for default###############
setwd("MACS_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS.run_i=data.motif.overlap.expected
    data.MACS.run_i$Run=i
    data.MACS.Allruns=rbind(data.MACS.Allruns,data.MACS.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                    Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                    Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                    Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                    Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                    Spatial.error.run_1=Spatial.error.run_1,
                                                    Spatial.error.run_2=Spatial.error.run_2,
                                                    Spatial.error.run_3=Spatial.error.run_3,
                                                    Spatial.error.run_4=Spatial.error.run_4,
                                                    Spatial.error.run_5=Spatial.error.run_5,
                                                    Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                              tot.motifs.run_2/tot.bs,
                                                                              tot.motifs.run_3/tot.bs,
                                                                              tot.motifs.run_4/tot.bs,
                                                                              tot.motifs.run_5/tot.bs)),
                                                    Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                              tot.motifs.run_2/tot.bs,
                                                                              tot.motifs.run_3/tot.bs,
                                                                              tot.motifs.run_4/tot.bs,
                                                                              tot.motifs.run_5/tot.bs)),
                                                    Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                tot.motifs.run_2/tot.bs,
                                                                                tot.motifs.run_3/tot.bs,
                                                                                tot.motifs.run_4/tot.bs,
                                                                                tot.motifs.run_5/tot.bs)),
                                                    Spatial.error.min=min(c(Spatial.error.run_1,
                                                                            Spatial.error.run_2,
                                                                            Spatial.error.run_3,
                                                                            Spatial.error.run_4,
                                                                            Spatial.error.run_5)),
                                                    Spatial.error.max=max(c(Spatial.error.run_1,
                                                                            Spatial.error.run_2,
                                                                            Spatial.error.run_3,
                                                                            Spatial.error.run_4,
                                                                            Spatial.error.run_5)),
                                                    Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                              Spatial.error.run_2,
                                                                              Spatial.error.run_3,
                                                                              Spatial.error.run_4,
                                                                              Spatial.error.run_5)),
                                                    tot.bs=x,Which="MACS")
    
    # merge:
    data.MACS.occurrence.spatial.error=rbind(data.MACS.occurrence.spatial.error,
                                             data.MACS.occurrence.spatial.error.x)
    
}
setwd("../../..")
getwd()
############ Load data MACS2 for default###############
setwd("MACS2_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS2.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS2.run_i=data.motif.overlap.expected
    data.MACS2.run_i$Run=i
    data.MACS2.Allruns=rbind(data.MACS2.Allruns,data.MACS2.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS2.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS2.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS2.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS2.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS2.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS2.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS2.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                     Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                     Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                     Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                     Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                     Spatial.error.run_1=Spatial.error.run_1,
                                                     Spatial.error.run_2=Spatial.error.run_2,
                                                     Spatial.error.run_3=Spatial.error.run_3,
                                                     Spatial.error.run_4=Spatial.error.run_4,
                                                     Spatial.error.run_5=Spatial.error.run_5,
                                                     Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                 tot.motifs.run_2/tot.bs,
                                                                                 tot.motifs.run_3/tot.bs,
                                                                                 tot.motifs.run_4/tot.bs,
                                                                                 tot.motifs.run_5/tot.bs)),
                                                     Spatial.error.min=min(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.max=max(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                               Spatial.error.run_2,
                                                                               Spatial.error.run_3,
                                                                               Spatial.error.run_4,
                                                                               Spatial.error.run_5)),
                                                     tot.bs=x,Which="MACS(2)")
    
    # merge:
    data.MACS2.occurrence.spatial.error=rbind(data.MACS2.occurrence.spatial.error,
                                              data.MACS2.occurrence.spatial.error.x)
    
}
setwd("../../..")
getwd()
############ Load data MACS3 for default###############
setwd("MACS3_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS3.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS3.run_i=data.motif.overlap.expected
    data.MACS3.run_i$Run=i
    data.MACS3.Allruns=rbind(data.MACS3.Allruns,data.MACS3.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS3.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS3.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS3.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS3.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS3.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS3.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS3.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                     Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                     Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                     Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                     Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                     Spatial.error.run_1=Spatial.error.run_1,
                                                     Spatial.error.run_2=Spatial.error.run_2,
                                                     Spatial.error.run_3=Spatial.error.run_3,
                                                     Spatial.error.run_4=Spatial.error.run_4,
                                                     Spatial.error.run_5=Spatial.error.run_5,
                                                     Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                 tot.motifs.run_2/tot.bs,
                                                                                 tot.motifs.run_3/tot.bs,
                                                                                 tot.motifs.run_4/tot.bs,
                                                                                 tot.motifs.run_5/tot.bs)),
                                                     Spatial.error.min=min(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.max=max(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                               Spatial.error.run_2,
                                                                               Spatial.error.run_3,
                                                                               Spatial.error.run_4,
                                                                               Spatial.error.run_5)),
                                                     tot.bs=x,Which="MACS(3)")
    
    # merge:
    data.MACS3.occurrence.spatial.error=rbind(data.MACS3.occurrence.spatial.error,
                                              data.MACS3.occurrence.spatial.error.x)
    
}
setwd("../../..")
getwd()
############ Load data MACS4 for default###############
setwd("MACS4_results")
setwd("Motifs")
setwd("Motifs.5000")
# initiate:
data.MACS4.Allruns=data.frame()
# load and merge data:
for(i in 1:5){
    setwd(paste("RUN_",i,"/Expected",sep=""))
    load("data.motif.overlap.expected.MACS")
    data.MACS4.run_i=data.motif.overlap.expected
    data.MACS4.run_i$Run=i
    data.MACS4.Allruns=rbind(data.MACS4.Allruns,data.MACS4.run_i)
    setwd("../..")
}
# now the data is loaded and merged
# create a global dataframe:
data.MACS4.occurrence.spatial.error=data.frame()
# fill in with values:
X.axis=seq(from=100,to=5000,by=100)#IF A2, A3
for(x in X.axis){
    tot.bs=x
    # Run 1:
    data.MACS.run_1=subset(data.MACS4.Allruns,Run==1)
    tot.motifs.run_1=nrow(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_1=mean(subset(data.MACS.run_1[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 2:
    data.MACS.run_2=subset(data.MACS4.Allruns,Run==2)
    tot.motifs.run_2=nrow(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_2=mean(subset(data.MACS.run_2[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 3:
    data.MACS.run_3=subset(data.MACS4.Allruns,Run==3)
    tot.motifs.run_3=nrow(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_3=mean(subset(data.MACS.run_3[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 4:
    data.MACS.run_4=subset(data.MACS4.Allruns,Run==4)
    tot.motifs.run_4=nrow(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_4=mean(subset(data.MACS.run_4[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # Run 5:
    data.MACS.run_5=subset(data.MACS4.Allruns,Run==5)
    tot.motifs.run_5=nrow(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=50))
    Spatial.error.run_5=mean(subset(data.MACS.run_5[1:x,],!is.na(Dist.to.expected.motif)&Dist.to.expected.motif<=100)$Dist.to.expected.motif)
    # -----create data frame:
    data.MACS4.occurrence.spatial.error.x=data.frame(Motif.occurence.run_1=tot.motifs.run_1/tot.bs,
                                                     Motif.occurence.run_2=tot.motifs.run_2/tot.bs,
                                                     Motif.occurence.run_3=tot.motifs.run_3/tot.bs,
                                                     Motif.occurence.run_4=tot.motifs.run_4/tot.bs,
                                                     Motif.occurence.run_5=tot.motifs.run_5/tot.bs,
                                                     Spatial.error.run_1=Spatial.error.run_1,
                                                     Spatial.error.run_2=Spatial.error.run_2,
                                                     Spatial.error.run_3=Spatial.error.run_3,
                                                     Spatial.error.run_4=Spatial.error.run_4,
                                                     Spatial.error.run_5=Spatial.error.run_5,
                                                     Motif.occurence.min=min(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.max=max(c(tot.motifs.run_1/tot.bs,
                                                                               tot.motifs.run_2/tot.bs,
                                                                               tot.motifs.run_3/tot.bs,
                                                                               tot.motifs.run_4/tot.bs,
                                                                               tot.motifs.run_5/tot.bs)),
                                                     Motif.occurence.mean=mean(c(tot.motifs.run_1/tot.bs,
                                                                                 tot.motifs.run_2/tot.bs,
                                                                                 tot.motifs.run_3/tot.bs,
                                                                                 tot.motifs.run_4/tot.bs,
                                                                                 tot.motifs.run_5/tot.bs)),
                                                     Spatial.error.min=min(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.max=max(c(Spatial.error.run_1,
                                                                             Spatial.error.run_2,
                                                                             Spatial.error.run_3,
                                                                             Spatial.error.run_4,
                                                                             Spatial.error.run_5)),
                                                     Spatial.error.mean=mean(c(Spatial.error.run_1,
                                                                               Spatial.error.run_2,
                                                                               Spatial.error.run_3,
                                                                               Spatial.error.run_4,
                                                                               Spatial.error.run_5)),
                                                     tot.bs=x,Which="MACS(4)")
    
    # merge:
    data.MACS4.occurrence.spatial.error=rbind(data.MACS4.occurrence.spatial.error,
                                              data.MACS4.occurrence.spatial.error.x)
    
}
setwd("../../..")
getwd()
############ Create the whole data ###############
Motif.occurence.allA3=rbind(data.MACPET.occurrence.spatial.error,data.MACS.occurrence.spatial.error,
                            data.MACS2.occurrence.spatial.error,
                            data.MACS3.occurrence.spatial.error,
                            data.MACS4.occurrence.spatial.error)
save("Motif.occurence.allA3",file=file.path(DIRsaveA1,"Motif.occurence.allA3"))
############ MO plot ###############

MOggplotA3=ggplot(Motif.occurence.allA3,aes(tot.bs,100*Motif.occurence.mean,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
# remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.y=element_blank(),
          axis.title.x=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0, 0),breaks=seq(from=0,to=100,by=1),limits=c(92,98))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.60),legend.title=element_blank(),
          legend.text=element_text(size=6))+guides(col = guide_legend(ncol=2))+
    scale_colour_manual(values=c("red","blue","green","black","purple"))#coloring factors


# MOggplotA3
save("MOggplotA3",file=file.path(DIRsaveA1,"MOggplotA3"))

############ SE plot ###############
SEggplotA3=ggplot(Motif.occurence.allA3,aes(tot.bs,Spatial.error.mean,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
# remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.y=element_blank())+xlab(" ")+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0, 0),breaks=seq(from=0,to=100,by=2),limits=c(10,18))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.95,.60),legend.title=element_blank(),
          legend.text=element_text(size=6))+guides(col = guide_legend(ncol=2))+
    scale_colour_manual(values=c("red","blue","green","black","purple"))#coloring factors

# SEggplotA3
save("SEggplotA3",file=file.path(DIRsaveA1,"SEggplotA3"))

############ FDR plot ###############
MACPET_run_1=subset(data.MACPET.Allruns,Run==1)
MACPET_run_1$FDR=MACPET_run_1$FDR*100/max(MACPET_run_1$FDR)
MACS_run_1=subset(data.MACS.Allruns,Run==1)
MACS_run_1$p.values.FDR=MACS_run_1$p.values.FDR*100/max(MACS_run_1$p.values.FDR)
MACS2_run_1=subset(data.MACS2.Allruns,Run==1)
MACS2_run_1$p.values.FDR=MACS2_run_1$p.values.FDR*100/max(MACS2_run_1$p.values.FDR)
MACS3_run_1=subset(data.MACS3.Allruns,Run==1)
MACS3_run_1$p.values.FDR=MACS3_run_1$p.values.FDR*100/max(MACS3_run_1$p.values.FDR)
MACS4_run_1=subset(data.MACS4.Allruns,Run==1)
MACS4_run_1$p.values.FDR=MACS4_run_1$p.values.FDR*100/max(MACS4_run_1$p.values.FDR)


data.FDRA3=data.frame(FDR=c(sort(MACPET_run_1$FDR),sort(MACS_run_1$p.values.FDR),
                            sort(MACS2_run_1$p.values.FDR),sort(MACS3_run_1$p.values.FDR),sort(MACS4_run_1$p.values.FDR)),
                      Which=c(rep("MACPET",nrow(MACPET_run_1)),rep("MACS",nrow(MACS_run_1)),
                              rep("MACS(2)",nrow(MACS2_run_1)),rep("MACS(3)",nrow(MACS3_run_1)),
                              rep("MACS(4)",nrow(MACS4_run_1))),
                      bs=c(1:nrow(MACPET_run_1),1:nrow(MACS_run_1),
                           1:nrow(MACS2_run_1),1:nrow(MACS3_run_1),1:nrow(MACS4_run_1)))

save("data.FDRA3",file=file.path(DIRsaveA1,"data.FDRA3"))

FDRggplotA3=ggplot(data.FDRA3,aes(bs,FDR,color=factor(Which),linetype=factor(Which)))+geom_line(size=0.3)+
    xlab("")+
    # remove background color, grid and add axes:
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8,family="Times"),
          axis.title=element_text(size=10,family="Times"),
          axis.title.y=element_blank())+
    ggplot2::scale_x_continuous(expand = c(0.01, 0))+ggplot2::scale_y_continuous(expand = c(0.005, 0),breaks=seq(from=0,to=100,by=10),limits=c(0,100))+
    # legend for MACPET and MACS
    theme(legend.justification=c(1,0), legend.position=c(.55,.30),legend.title=element_blank(),
          legend.text=element_text(size=6))+
    scale_colour_manual(values=c("red","blue","green","black","purple"))#coloring factors

# FDRggplotA3
save("FDRggplotA3",file=file.path(DIRsaveA1,"FDRggplotA3"))