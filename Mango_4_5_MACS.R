###################################################################################################################
#################################### Script for MACS-MANGO in R ################################################
###################################################################################################################
# Important: One has to specify a save directory and a data directory, the current module only gives information about how the models were ran.
library(MACPET)
#---------------------
#Specify the folders and the data:
#---------------------
DATA.type="ChIA_PET_Data"
Experiment="ENCSR000FDG"#experiment ID, other experiments to choose: ENCSR000CAD, ENCSR000CAC, ENCSR000FDD, ENCSR000BZZ, ENCSR000BZY
Repetition="REP1"
DIRsave=file.path("YOURPATH",DATA.type,Experiment,Repetition,"InteractionsResults_MANGO_MACPETInitials/MACS_N")# the data directory
if(!dir.exists(DIRsave)) dir.create(DIRsave)
#---------------------
#Define where the peaks are:
#---------------------
Load.MACS.Peaks=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACS_results/MACS.peaks.out_peaks.update")
load(Load.MACS.Peaks)
# keep significant:--------------------NOTE, here we have teste 0.05,0.01 and N=total significant from MACPET as requested. The total significant from MACPET can be found by uploading the psfit file of MACPET of the corresponding dataset and count nrows(file) on the 0.05 threshold.
MACS.peaks.out_peaks.update=subset(MACS.peaks.out_peaks.update,p.values.FDR<0.05)
# make GRanges:
#make GRanges:
MACS_peaks.narrowPeak=GenomicRanges::GRanges(seqnames=MACS.peaks.out_peaks.update$chr,
                           ranges=IRanges::IRanges(
                               start=round(MACS.peaks.out_peaks.update$start),
                               end=round(MACS.peaks.out_peaks.update$end)))
#----sort by chrom and ranges:
MACS_peaks.narrowPeak=GenomeInfoDb::sortSeqlevels(MACS_peaks.narrowPeak)
MACS_peaks.narrowPeak=S4Vectors::sort(MACS_peaks.narrowPeak)
#-----add names:
MACS_peaks.narrowPeak$name=paste("peak_",1:length(MACS_peaks.narrowPeak),sep="")
# save:
connection=file.path(DIRsave,"MACS_peaks.narrowPeak")#since MANGO needs this like that
rtracklayer::export(MACS_peaks.narrowPeak,con=connection,format="BED")
#---------------------
#-----All stages parameters
#---------------------
MangoPath="YOURPATH/mango/mango.R"#the R file for the mango, in the external drive
stages="4:5"#IF i use stage 4 here it will extend and merge the peaks, try not to do that too.
prefix="MACS"
chromexclude="chrM,chrY"#to exclude always
chrominclude="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrM,chrY"
#macs2 path:
macs2path="/usr/local/bin/macs2"
#the bowtie genome reference
bowtieref= "YOURPATH/Bowtie/hg19.ebwt/hg19"#in the external drive
#the bedtools genome reference
bedtoolsgenome="YOURPATH/bedtools/bedtools2/genomes/human.hg19.genome"#in the external drive
#bowtie path:
bowtiepath="YOURPATH/Bowtie/bowtie-1.2/bowtie"#in the external
#bedtools path:
bedtoolspath="YOURPATH/bedtools/bedtools2/bin/bedtools"#in the external
#outdir:
outdir=DIRsave#where to save everything
#---------------------
#---------- STAGE 4 PARAMETERS (peakslop is defined inside the loop.)
#---------------------
#----Only if this is run by mango, else provide your own peaks
MACS_qvalue="0.05"#pvalue cutoff for peak calling in MACS2. default = 0.05#NOT USED BUT MANGO NEEDS IT AS INPUT ANYWAY
peakinput=connection#Name of user supplied peaks file. If NULL Mango will use peaks determined from MACS2 analysis. default = NULL
MACS_shiftsize="NULL"#let MACS2 decide#NOT USED BUT MANGO NEEDS IT AS INPUT ANYWAY
blacklist="NULL"#NULL, removed
gsize="hs"#mappable genome size or effective genome size for MACS2.default = 'hs'#NOT USED BUT MANGO NEEDS IT AS INPUT ANYWAY
#---------------------
#---------- STAGE 5 PARAMETERS
#---------------------
# Load the self ligated cut:
DIRMACPETdata=file.path("YOURPATH",DATA.type,Experiment,Repetition,"MACPET_results/S3_results/MACPET_psfitData")
load(DIRMACPETdata)
selfcut=as.character(metadata(MACPET_psfitData)$MaxSize)#My self-ligation cutoff! Because you need to use the same self-ligation cut-off as MACPET for the corresponding dataset. To ensure that the intra-chromosomal are the same as the ones used in MACPET
# #note to get the same significant as MACPET for MACS, that is N:
# nrow(subset(metadata(MACPET_psfitData)$Peaks.Info,FDR<0.05))
# MACS.peaks.out_peaks.update=MACS.peaks.out_peaks.update[with(MACS.peaks.out_peaks.update,order(p.values.FDR,decreasing=F)),]
# nrow(MACS.peaks.out_peaks.update)
# MACS.peaks.out_peaks.update=MACS.peaks.out_peaks.update[1:nrow(subset(metadata(MACPET_psfitData)$Peaks.Info,FDR<0.05)),]
# nrow(MACS.peaks.out_peaks.update)
# #end note
distcutrangemin="1000" #When Mango determines the self-ligation cutoff this is the minimum distance it will consider. Changing this setting is not recommended. default = 1000#NOT USED BUT MANGO NEEDS IT AS INPUT ANYWAY
distcutrangemax="100000"#When Mango determines the self-ligation cutoff this is the maximum distance it will consider. Changing this setting is not recommended. default = 100000#NOT USED BUT MANGO NEEDS IT AS INPUT ANYWAY
biascut="0.05" #Mango exlcudes very short distance PETS since they tend to arise from self-ligation of a single DNA framgent as opposed to interligation of two interacting fragments. To determine this distnce cutoff Mango determines the fraction of PETs at each distance that come from self-ligation and sets the cutoff at the point where the fraction is less than or equal to BIASCUT. default = 0.05
FDR="0.05"#FDR cutoff for significant interactions. default = 0.05
numofbins="50"#number of bins to use for binomial p-value calculations. default = 50
corrMethod="BH" #Method to use for correction of mulitply hypothesis testing. See (http://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html) for more details. default = BH
maxinteractingdist="1000000"#The maximum disance (in basepairs) considered for interaction. Optimum sensitivity is generally acheived at values of 1000000-2000000. default = 1000000
extendreads="120"#how many bp to extend reads towards peak. default = 120
minPETS="2"#The minimum number of PETs required for an interaction (applied after FDR filtering). default = 2
reportallpairs="FALSE"#Should all pairs be reported or just significant pairs (TRUE or FALSE). default = FALSE
#---------------------
#----------Create the command
#---------------------
#first you need to buildt the tagalign files: Because MANGO cannot actually run from stage 4, although they say they do..
outname = file.path( outdir ,prefix)
bedpefilesortrmdup = paste(outname ,".rmdup.bedpe",sep="")#This file is the intra-chromosomal from MACPET, converted to bedpe
tagAlignfile       = paste(outname,".tagAlign",sep="")
mango::buildTagAlign(bedpefilesortrmdup ,tagAlignfile )
#---------------------
#----------run a loop command for all the sizes, 
#---------------------
Peak_merge_seq=seq(0,1000,100)
for(i in Peak_merge_seq){
    peakslop=as.character(i)#the peak extention., the other stage 4 parameter
    # create the commant:
    Command=paste("Rscript",MangoPath,"--stages",stages,
                  "--chrominclude",chrominclude,"--chromexclude",chromexclude,
                  "--prefix",prefix,"--outdir",outdir,"--bowtieref",bowtieref,
                  "--bedtoolsgenome",bedtoolsgenome,"--bowtiepath",bowtiepath,
                  "--bedtoolspath",bedtoolspath,"--macs2path",macs2path,
                  "--MACS_qvalue",MACS_qvalue,"--peakslop",peakslop,"--peakinput",peakinput,
                  "--MACS_shiftsize",MACS_shiftsize,"--blacklist",blacklist,
                  "--gsize",gsize,"--distcutrangemin",distcutrangemin,
                  "--distcutrangemax",distcutrangemax,"--biascut",biascut,
                  "--FDR",FDR,"--numofbins",numofbins,"--corrMethod",corrMethod,
                  "--maxinteractingdist",maxinteractingdist,"--extendreads",extendreads,
                  "--minPETS",minPETS,"--reportallpairs",reportallpairs,"--selfcut",selfcut)
    # run:
    system(Command)
    # move files to the given directory:
    From=outdir
    To=file.path(outdir,paste("Interactions_",i,sep=""))
    if(!dir.exists(To)) dir.create(To)
    # 1:
    file.rename(from=file.path(From,"MACS.models.pdf"),to=file.path(To,"MACS.models.pdf"))
    # 2:
    file.rename(from=file.path(From,"MACS_peaks.slopPeak"),to=file.path(To,"MACS_peaks.slopPeak"))
    # 3:
    file.rename(from=file.path(From,"MACS_peaks.slopPeak.depth"),to=file.path(To,"MACS_peaks.slopPeak.depth"))
    # 4:
    file.rename(from=file.path(From,"MACS.depth_combo_model.1.text"),to=file.path(To,"MACS.depth_combo_model.1.text"))
    # 5:
    file.rename(from=file.path(From,"MACS.depth_combo_model.2.text"),to=file.path(To,"MACS.depth_combo_model.2.text"))
    # 6:
    file.rename(from=file.path(From,"MACS.depth_IAB_model.1.text"),to=file.path(To,"MACS.depth_IAB_model.1.text"))
    # 7:
    file.rename(from=file.path(From,"MACS.depth_IAB_model.2.text"),to=file.path(To,"MACS.depth_IAB_model.2.text"))
    # 8:
    file.rename(from=file.path(From,"MACS.distance_combo_model.1.text"),to=file.path(To,"MACS.distance_combo_model.1.text"))
    # 9:
    file.rename(from=file.path(From,"MACS.distance_combo_model.2.text"),to=file.path(To,"MACS.distance_combo_model.2.text"))
    # 10:
    file.rename(from=file.path(From,"MACS.distance_IAB_model.1.text"),to=file.path(To,"MACS.distance_IAB_model.1.text"))
    # 11:
    file.rename(from=file.path(From,"MACS.distance_IAB_model.2.text"),to=file.path(To,"MACS.distance_IAB_model.2.text"))
    # 12:
    file.rename(from=file.path(From,"MACS.interactions.fdr.mango"),to=file.path(To,"MACS.interactions.fdr.mango"))
    # 13:
    file.rename(from=file.path(From,"MACS.mango.log"),to=file.path(To,"MACS.mango.log"))
}


