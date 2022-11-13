
#-----------------------------------------------------------------
# simulate the source population with non-overlapping generations
#-----------------------------------------------------------------
setwd("~/Documents/orca/analyses_31May2022/Projections")

# designate historical effective population size
pastNe <- round(read.table("SRKWHistoricalNe",header=TRUE)[1:150,2])
pastNe <- c(pastNe,rep(pastNe[length(pastNe)],849))
N <- c(rev(pastNe),c(120,120))          
Gens <- 1000        # total number of generations to run the source population simulation for
Shape <- 0.2        # shape parameter for the gamma distribution of selection coefficients
Scale=0.2           # scale parameter for the gamma distribution of selection coefficients
PropLethal <- 0.01  # proportion of mutations to make lethal when homozygous
beta <- 50          # exponential rate of decreasing dominance coefficient (h) with increasing size of the selection coefficient (s). 
                    # The dominance model is from Deng & Lynch (1996, Genetics 144 349-360).
DelMu <- 9e-7       # deleterious mutation rate per base pair in the simulated genome. Note that the simulated genome is quite small physically (genSize parameter below), but this size is completely arbitrary. 
library(quantPop)
logQuant_mut_delMut_ceiling(burnin=1,
                            gens=1001,
                            genSize=10000000,
                            chroms=20,
                            phen0=0,
                            phenOpt=rep(0,1001),
                            c=6,
                            mu=0,
                            quantDom = FALSE,
                            muOff=1001,
                            minSize=-0.5,
                            maxSize=0.5,
                            N=N,
                            Ve=4,
                            hardSelecGen=1002,
                            K=N,
                            lambda=2.5,
                            f=4,
                            delMutation=TRUE,
                            delMu=DelMu,
                            propLethal=PropLethal,
                            delMuGammaShape=Shape,
                            delMuGammaScale=Scale,
                            neutMutation=FALSE,
                            neutMu=0,
                            Beta=beta,
                            importGenos=FALSE,
                            importGenoIndivs=NULL,
                            mutationTag=NULL,
                            genoMat1Name=NULL,
                            genoMat2Name=NULL,
                            locusInfoName=NULL,
                            genRescue=FALSE,
                            rescueInGenos1=NULL,
                            rescueInGenos2=NULL,
                            rescueLocusInfo=NULL,
                            rescueN=NULL,
                            rescueGens=NULL,
                            pedigree=FALSE,
                            multModel = TRUE,
                            saveFreqs = FALSE)



# write the input
write.table(genoMat1,file=paste("inputKillerGenos1_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(genoMat2,file=paste("inputKillerGenos2_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(polyInfo,file=paste("inputKillerLocusInfo_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
