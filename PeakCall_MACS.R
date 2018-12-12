####################################################################################################################
####################################                 MACS                 ##########################################
####################################################################################################################
# Module for running MACS on the different parameter setings
# The same module was used for the default settings too
# Important: One has to specify a save directory and a data directory, the current module only gives information about how the models were ran.
library(MACPET)
library(InteractionSet)
####################################################################################################################
################## Those are MACS at the default parameters bw=300, slocal##########################################
####################################################################################################################
#------------------------------------------------------------------
#---------------------------------Running:
#------------------------------------------------------------------
#-------
#-----Data preperation:
#------
Experiment="ENCSR000FDG"#experiment ID, other experiments to choose: ENCSR000CAD, ENCSR000CAC, ENCSR000FDD, ENCSR000BZZ, ENCSR000BZY
Repetition="REP1"#repetition of the experiment
AnalysisDir=file.path("YOURPATH",Experiment,Repetition)
load(paste(AnalysisDir,"/MACPET_results/S3_results/MACPET_pselfData",sep=""))#load self ligated PETs from MACPET
pselfData=MACPET_pselfData
# Create a save dir for the MACS output:
DIRMACS=file.path(AnalysisDir,"MACS_results")
if(!dir.exists(DIRMACS)) dir.create(DIRMACS)
MACSdir="YOURPATH/MACS-1.4.2/build/scripts-2.7/macs14"#where MACS is installed
#-------
#-----Data preperation:
#------
#-------------Keep 5 ends as MACS says:
Anchor1=InteractionSet::anchors(pselfData,type="first")
Anchor1=methods::as(Anchor1,"GRanges")
# invert strands since GRanges has inverted them (or else MACS will find no peaks):
SP=which(strand(Anchor1)=="+")
SM=which(strand(Anchor1)=="-")
strand(Anchor1[SP])="-"
strand(Anchor1[SM])="+"
#save for MACS to get as input:
con=paste(DIRMACS,"/Anchor1",sep="")
rtracklayer::export(Anchor1,con=con,format="BED")
#--------------------------------------------
#-----Command creation: RUN ONE OF THE FOLLOWING MODELS for each parameter setting for MACS.
# Note that for the other analysis with the different thresholds you dont have to run MACS again, just take a lower threshold.
#--------------------------------------------

################## Those are MACS at the default parameters bw=300, slocal  (1)
SaveFile=paste(DIRMACS,"/MACS.peaks.out",sep="")
Command=paste(MACSdir, "-p 0.05 -f AUTO -t", con, "-n",SaveFile)
system(Command)#run it

###############  Those are MACS at bw=600, slocal (2)
SaveFile=paste(DIRMACS,"/MACS.peaks.out",sep="")
Command=paste(MACSdir, "-p 0.05 --bw 600 -f AUTO -t", con, "-n",SaveFile,"--slocal 1000")
system(Command)#run it

############## Those are MACS at bw=1000, llocal (3)
SaveFile=paste(DIRMACS,"/MACS.peaks.out",sep="")
Command=paste(MACSdir, "-p 0.05 --bw 1000 -f AUTO -t", con, "-n",SaveFile,"--llocal 10000")
system(Command)#run it


############### Those are MACS at bw=300, llocal (4)
SaveFile=paste(DIRMACS,"/MACS.peaks.out",sep="")
Command=paste(MACSdir, "-p 0.05 --bw 300 -f AUTO -t", con, "-n",SaveFile,"--llocal 10000")

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------- After running MACS for eithe model------------------------------------
#-------------------------------------------------------------------------------------------------------------------
setwd(DIRMACS)
#----------------
#----------Fix results in a better form:
#----------------
MACS.peaks.out_peaks.update=utils::read.table("MACS.peaks.out_peaks.xls",stringsAsFactors=FALSE)
ColNames=as.character(MACS.peaks.out_peaks.update[1,])
MACS.peaks.out_peaks.update=MACS.peaks.out_peaks.update[-1,]
MACS.peaks.out_peaks.update$V1=as.character(MACS.peaks.out_peaks.update$V1)#chrom
MACS.peaks.out_peaks.update$V2=as.numeric(MACS.peaks.out_peaks.update$V2)#start
MACS.peaks.out_peaks.update$V3=as.numeric(MACS.peaks.out_peaks.update$V3)#end
MACS.peaks.out_peaks.update$V4=as.numeric(MACS.peaks.out_peaks.update$V4)#size
MACS.peaks.out_peaks.update$V5=as.numeric(MACS.peaks.out_peaks.update$V5)#summit from +start
MACS.peaks.out_peaks.update$V6=as.numeric(MACS.peaks.out_peaks.update$V6)#tags
MACS.peaks.out_peaks.update$V7=as.numeric(MACS.peaks.out_peaks.update$V7)#p-value -log
MACS.peaks.out_peaks.update$V8=as.numeric(MACS.peaks.out_peaks.update$V8)#fold
colnames(MACS.peaks.out_peaks.update)=c("chr","start","end","length","summit","tags","p.values.log","fold_enrichment")
#transform p and q values as documented in MACS:
MACS.peaks.out_peaks.update$p.values=exp(MACS.peaks.out_peaks.update$p.values.log*log(10)/(-10))
MACS.peaks.out_peaks.update$p.values.FDR=p.adjust(p=MACS.peaks.out_peaks.update$p.values,method = "BH")
save(MACS.peaks.out_peaks.update,file="MACS.peaks.out_peaks.update")
