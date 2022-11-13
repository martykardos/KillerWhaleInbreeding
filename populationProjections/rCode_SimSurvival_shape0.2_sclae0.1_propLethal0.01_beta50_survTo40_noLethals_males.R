# cleanup
rm(list=ls())
gc()

setwd("~/Documents/orca/analyses_31May2022/Projections")

# survival model is truncated for younger animals.
# input mutation parameters
Shape <- 0.2
Scale=0.1
PropLethal <- 0.01
beta <- 50
SimSex <- "male"
B <- 3.74 # haploid lethal equivalents for male survival to 40 years (the empirical observation)
reps <- 50

library(quantPop)
repIndex <- 1    # which set of simulations are these
NYEARS =  100 # number of years to project
KVec <- rep(200,NYEARS) # vector of carrying capacities
nSims = 500 # number of simulations
fullRecessive <- FALSE                     # run the simulations using fully recessive deleterious alleles?
simInbDep  <- TRUE      # do you want to include inbreeding depression?

#--------------------------------------------------------------------------
# read the input genomes from the source population
#--------------------------------------------------------------------------
inGenos1 <- read.table(paste("inputKillerGenos1_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,sep=""),header=FALSE)
inGenos2 <- read.table(paste("inputKillerGenos2_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,sep=""),header=FALSE)
inLocusInfo <- read.table(paste("inputKillerLocusInfo_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,sep=""),header=FALSE)
colnames(inGenos1) <- inLocusInfo[,2]
colnames(inGenos2) <- inLocusInfo[,2]

# remove lethals
PropLethal <- 0.0    # switch to proportion of lethals = 0 for the simulations below
lethals <- which(inLocusInfo[,4] == "1")
inGenos1 <- inGenos1[,-lethals]
inGenos2 <- inGenos2[,-lethals]
inLocusInfo <- inLocusInfo[-lethals,]

# calculate the inbreeding load
sampB <- NULL
p <- colSums(rbind(inGenos1,inGenos2))/(2*nrow(inGenos1))
if(fullRecessive == FALSE) hCoef <- 0.5*(exp(-beta*as.numeric(inLocusInfo[,4])))
if(fullRecessive == TRUE) hCoef <- rep(0,nrow(inLocusInfo))
singleLocusB <- p*as.numeric(inLocusInfo[,4]) - (p^2)*as.numeric(inLocusInfo[,4]) - 2*(p*(1-p)*as.numeric(inLocusInfo[,4])*hCoef)
sampB<- sum(p*as.numeric(inLocusInfo[,4])) - sum((p^2)*as.numeric(inLocusInfo[,4])) - 2*sum(p*(1-p)*as.numeric(inLocusInfo[,4])*hCoef)
print(sampB)

# remove loci that are fixed
remFixed <- which(p == 1)
inGenos1 <- inGenos1[,-remFixed]
inGenos2 <- inGenos2[,-remFixed]
inLocusInfo <- inLocusInfo[-remFixed,]

# libraries
library(lattice)
library(quantPop)
library(kinship2)
library(scales)
#preliminaries

survivModel <- readRDS("fit.rds")
# load in demographic data and parameters
load("demog_models_projection.Rdata")
# scenario used affects which rates are used for projections
years_to_sample = 2000

# These were the sex ratios we'd used during the workshops -- this is
# a random variable, rather than fixed, so the uncertainty is propagated
pFemale = 50/100 # this is the number of female births/total up through 2018
sdFemale = sqrt(pFemale*(1-pFemale)/100)

currentPop$froh = rep(-1.2622,nrow(currentPop)) # in reality, raw values were standardized using mean = 0.03234788, sd = 0.02563815


#store outputs
loadMat<- matrix(NA,nrow=nSims,ncol=NYEARS+1)               # total genetic load
meanFMat <- matrix(NA,nrow=nSims,ncol=NYEARS+1)
meanHMat <- matrix(NA,nrow=nSims,ncol=NYEARS+1)
popSize = array(0, dim=c(nSims, NYEARS))
BMat <- matrix(NA,nrow=nSims,ncol=NYEARS+1)

# to make this more general, I called "populationSegment" = pod, and "breedingGroup" = matriline
initPopCurrent = currentPop[,c("animal", "pod", "matriline", "sexF1M2", "age","froh")] # animal, pod, matriline, sex, age
initPopCurrent = dplyr::rename(initPopCurrent, sex=sexF1M2)
initPopCurrent$animal = as.character(initPopCurrent$animal)
initPopCurrent$matriline = as.numeric(as.factor(initPopCurrent$matriline)) # , -29 because of NRs
initPopCurrent$dad = NA # keep track of dads

numWhale = dim(initPopCurrent)[1]
initPopCurrent$hasCalf = 0
initPopCurrent$hasCalf[which(as.character(currentPop$animal) %in% as.character(currentPop$mom[which(currentPop$age==1)]))] = 1
initPopCurrent$mom = as.character(currentPop$mom)

initPopCurrent$sex = as.numeric(initPopCurrent$sex)
# assign sexes to unsexed individuals [0s]
initPopCurrent$sex[which(initPopCurrent$sex==0)] = ifelse(runif(length(which(initPopCurrent$sex==0))) < rnorm(length(which(initPopCurrent$sex==0)),pFemale, sdFemale), 1, 2)
newID = 9999

#---------------------------------------------------
# read in the Southern Resident pedigree
#---------------------------------------------------
SRPed <- read.table("pedigreeForSimulations_30Nov2021",header=TRUE)

# add new individuals that are in the population database, with only maternal information available, to the pedigree
onlyMomIDs <- initPopCurrent[which(initPopCurrent[,1] %in% SRPed[,1] == FALSE),1]

#add moms that do not have a record in the pedigree (i.e., a row committed to them)
missingMoms <- unique(initPopCurrent[which(initPopCurrent[,9] %in% SRPed[,1] == FALSE),9])
missingMoms <- missingMoms[-which(is.na(missingMoms))]
missingMoms <- missingMoms[-which(missingMoms == "L091")]

# updated pedigree
newPed <- matrix(NA,nrow=length(onlyMomIDs) + length(missingMoms),ncol=ncol(SRPed))
newPed[,1] <- c(onlyMomIDs,missingMoms)
newPed[1:length(onlyMomIDs),6] <- initPopCurrent[match(onlyMomIDs,initPopCurrent[,1]),9]   # moms of individuals with only mom ids


newPed[1:length(onlyMomIDs),2] <- c(rep("2020",4),"1985",rep("2020",6))                    # say that they are born in 2020 to ensure the moms are alive when simulating their genomes
newPed[1:length(onlyMomIDs),4] <- initPopCurrent[match(onlyMomIDs,initPopCurrent[,1]),4]   # sex of individuals with only mom IDs available
newPed[(length(onlyMomIDs)+1):(nrow(newPed)),4] <- 1                                       # moms are all female


newPed <- as.data.frame(newPed)
colnames(newPed) <- colnames(SRPed)
SRPed <- rbind(SRPed,newPed)

# add founders to accomodate mating to produce individuals with only one known parent. Assume the other parent is randomly drawn from the source population.
addFounderDads <- sum(is.na(SRPed$dad) == TRUE & is.na(SRPed$mom) == FALSE)
addFounderMoms <- sum(is.na(SRPed$mom) == TRUE & is.na(SRPed$dad) == FALSE)
newDadFoundPedInfo <- matrix(NA,nrow=addFounderDads,ncol=ncol(SRPed))
newDadFoundPedInfo[,1] <- c(paste("M000",6:9,sep=""),paste("M00",10:23,sep=""))
newDadFoundPedInfo[,4] <- rep(2,addFounderDads)

newMomFoundPedInfo <- matrix(NA,nrow=addFounderMoms,ncol=ncol(SRPed))
newMomFoundPedInfo[,1] <- c(paste("F000",4:9,sep=""),paste("F00",10:14,sep=""))
newMomFoundPedInfo[,4] <- rep(1,addFounderMoms)

colnames(newDadFoundPedInfo) <- colnames(SRPed)
colnames(newMomFoundPedInfo) <- colnames(SRPed)

# add new dad and mom ids to the father and mother column of the pedigree
SRPed$dad[which(is.na(SRPed$dad) == TRUE & is.na(SRPed$mom) == FALSE)] <- c(paste("M000",6:9,sep=""),paste("M00",10:23,sep=""))
SRPed$mom[which(is.na(SRPed$mom) == TRUE & is.na(SRPed$dad) == FALSE)] <- c(paste("F000",4:9,sep=""),paste("F00",10:14,sep=""))

SRPed <- rbind(SRPed,newDadFoundPedInfo,newMomFoundPedInfo)

# save output
outSurvProbs <- matrix(NA,nrow=reps,ncol=74)
outSurvProbs2 <- matrix(NA,nrow=reps,ncol=74)
outFh <- matrix(NA,nrow=reps,ncol=74)
outHet <- matrix(NA,nrow=reps,ncol=74)
realizedB <- rep(NA,reps)
for(y in 1:reps){
    #---------------------------------------------------------------------------------------------------
    # randomly assign each pedigree founder a genome from the source population
    #---------------------------------------------------------------------------------------------------
    founders <- SRPed[is.na(SRPed$dad) & is.na(SRPed$mom),]
    nonFounders <- SRPed[which((is.na(SRPed$dad) & is.na(SRPed$mom)) == FALSE),]
    sampGens <- sort(sample(1:nrow(inGenos1),nrow(founders),replace=FALSE))    # which genomes to assign to the pedigree founders
    #delLoci <- 1:nrow(inGenos1)                                                 # which loci have deleterious effects on viability
    
    newLocusInfo <- as.matrix(inLocusInfo)
    newLocusInfo[,1] <- rep(0,nrow(newLocusInfo))
    
    
    genos1 <- inGenos1[sampGens,]                               # first haploid genome copies
    genos2 <- inGenos2[sampGens,]                               # second haploid genome copies
    colnames(genos1) <- as.character(newLocusInfo[,2])
    colnames(genos2) <- as.character(newLocusInfo[,2])
    # add founder genomes that will mate with a known parent to produce offspring with only one known parent
    genosIDs <- founders[,1]                              # link IDs to genomes
    
    
    
    
    #--------------------------------------------------------------------------------------------------------------------------
    # prune the unnecessary loci that are present in the founders so that the initial inbreeding load for survival to 40 years is B ~ 8
    # This is to make simulating the pedigree faster/
    #--------------------------------------------------------------------------------------------------------------------------
    
    B0 <- 1.1*B
    sampB <- NULL
    p <- colSums(rbind(genos1,genos2))/(2*nrow(genos1))
    if(fullRecessive == FALSE) hCoef <- 0.5*(exp(-beta*as.numeric(newLocusInfo[,4])))
    if(fullRecessive == TRUE) hCoef <- rep(0,nrow(newLocusInfo))
    singleLocusB <- p*as.numeric(newLocusInfo[,4]) - (p^2)*as.numeric(newLocusInfo[,4]) - 2*(p*(1-p)*as.numeric(newLocusInfo[,4])*hCoef)
    sampB<- sum(p*as.numeric(newLocusInfo[,4])) - sum((p^2)*as.numeric(newLocusInfo[,4])) - 2*sum(p*(1-p)*as.numeric(newLocusInfo[,4])*hCoef)
    
    
    singleLocusB2 <- singleLocusB
    newLocusInfo2 <- newLocusInfo
    bVec <- NULL
    nLocusVec <- NULL
    if(sampB >= B0){
      testThisB <- FALSE
      while(testThisB == FALSE){
        if(sum(singleLocusB2) >= B0){
          candLocus <- sample(1:nrow(newLocusInfo2),1,replace=FALSE)
          newLocusInfo2 <- newLocusInfo2[-candLocus,]
          singleLocusB2 <- singleLocusB2[-candLocus]
        }
        if(sum(singleLocusB2) <= B0){
          testThisB <- TRUE
        }
        bVec <- c(bVec,sum(singleLocusB2))
        nLocusVec <- c(nLocusVec,length(singleLocusB2))
        #print(paste("pruned load = ",sum(singleLocusB2) ))
      }
      keepLoci <- match(newLocusInfo2[,2],newLocusInfo[,2])
      
      genos1 <- genos1[,keepLoci]
      genos2 <- genos2[,keepLoci]
      newLocusInfo <- newLocusInfo[keepLoci,]
    }
    
    hets <- rowSums((genos1 != genos2))/ncol(genos1)
    if(fullRecessive == FALSE) getIndivLoads(genomeInfo=newLocusInfo,inGenomes1=genos1,inGenomes2=genos2,expRate=beta,fullRecessive=FALSE,multModel=TRUE)
    if(fullRecessive == TRUE) getIndivLoads(genomeInfo=newLocusInfo,inGenomes1=genos1,inGenomes2=genos2,expRate=beta,fullRecessive=TRUE,multModel=TRUE)
    
    # set the reference genetic load as that expected for individual with the highest heterozygosity (assumed to be nob-inbred). This works because we are modelling inbreeding depression within the population. Drift load is assumed fixed and reflected in the fitness of all individuals in the population.
    Fh <- (max(hets)-hets)/max(hets) 
    mod <- lm(indivLoads~Fh)
    refLoad <- mod[[1]][1]

    
    #---------------------------------------------------
    # simulated mating in the known pedigree
    #---------------------------------------------------
    #offspring <- SRPed[which(SRPed[,1] %in% genosIDs == FALSE),1]
    offYears <- sort(unique(nonFounders[,2]))
    for(j in offYears){
      theseOffIDs <-  nonFounders[nonFounders[,2] == j,1]
      theseDads <- SRPed$dad[match(theseOffIDs,SRPed$ID)]
      theseMoms <- SRPed$mom[match(theseOffIDs,SRPed$ID)]
      theseSexes <- SRPed$M1F2[match(theseOffIDs,SRPed$ID)]
      mutRate <- 5e-8      # set the per bp mutation rate so you have one deleterious mutation per average offspring
      
      for(kl in 1:length(theseOffIDs)){
        survived <- 0 # indicator for whether the offspring survives
        while(survived == 0){
          # make a candidate offspring genome. It has to survive to be included in the pedigree
          mating_delMutation(delMu = mutRate,
                             genSize=10000000,
                             chroms=20,
                             mothers= theseMoms[kl],
                             fathers= theseDads[kl],
                             genoIDVec = genosIDs,
                             propLethal=PropLethal,
                             delMuGammaShape=Scale,
                             delMuGammaScale=Shape,
                             Beta=beta,
                             importGenos=TRUE,
                             genoMat1Name=genos1,
                             genoMat2Name=genos2,
                             inLocusInfo=newLocusInfo,
                             year=1)
          #------------------------------------------
          # see if the candidate offspring survived
          #------------------------------------------
          if(fullRecessive == FALSE) getIndivLoads(genomeInfo=polyInfo,inGenomes1=genoMat1,inGenomes2=genoMat2,expRate=beta,fullRecessive=FALSE,multModel=TRUE)   # the inputs here are the genotypes of the simulated offspring arising from the call to mating_delMutation() just above here
          if(fullRecessive == TRUE) getIndivLoads(genomeInfo=polyInfo,inGenomes1=genoMat1,inGenomes2=genoMat2,expRate=beta,fullRecessive=TRUE,multModel=TRUE)
      
          if(indivLoads != 1) scaleLoads <- indivLoads - refLoad      # scale loads relative to the individual with the highest heterozygosity
          if(indivLoads == 1) scaleLoads <- 1      # allow lethal homozygous genotypes to be lethal

          # calculate reference survival to reproductive age probability
          survProbRaw <- prod(c(0.5,rep(0.9934952,10)))    # 0.9956314 is the yearly survival probability of a minimally inbred adult female. Note that prod(rep(0.9956314,39)) == 0.8430335, which is nearly identical to the inferred probability of a non-inbred female surviving to 40 (see supplementary demographic figures)

          # incorporate inbreeding depression on survival to reproductive age
          survProbLoad <- survProbRaw - scaleLoads*survProbRaw
          survTest <- runif(1,min=0,max=1) <= survProbLoad      # does the simulated individual survive to 10 years?
          if(survTest == TRUE)survived <- 1                     # continue if the offspring survived
        }
        
        if(ncol(genoMat1) == ncol(genos1)){
          genos1 <- rbind(genos1,genoMat1)          # update genotype matrics and locus information with new individuals
          genos2 <- rbind(genos2,genoMat2)
        }
        
        if(ncol(genoMat1) != ncol(genos1)){
          mutLocs <- as.numeric(polyInfo[which(polyInfo[,6] %in% newLocusInfo[,6] == FALSE),6])
          newLocusLocs <- as.numeric(newLocusInfo[,6])
          for(k in 1:length(mutLocs)){
            if(mutLocs[k] < min(newLocusLocs)){
              genos1 <- cbind(rep(FALSE,nrow(genos1)),genos1)
              genos2 <- cbind(rep(FALSE,nrow(genos2)),genos2)
              newLocusLocs <- c(mutLocs[k],newLocusLocs)
            }
            if(mutLocs[k] > min(as.numeric(newLocusInfo[,6])) & mutLocs[k] < max(as.numeric(newLocusInfo[,6]))){
              genos1 <- cbind(genos1[,which(newLocusLocs < mutLocs[k])],rep(FALSE,nrow(genos1)),genos1[,which(newLocusLocs > mutLocs[k])] )
              genos2 <- cbind(genos2[,which(newLocusLocs < mutLocs[k])],rep(FALSE,nrow(genos2)),genos2[,which(newLocusLocs > mutLocs[k])] )
              newLocusLocs <- c(newLocusLocs[which(newLocusLocs < mutLocs[k])],mutLocs[k],newLocusLocs[which(newLocusLocs > mutLocs[k])])
            }
            if(mutLocs[k] > max(as.numeric(newLocusInfo[,6]))){
              genos1 <- cbind(genos1,rep(FALSE,nrow(genos1)))
              genos2 <- cbind(genos2,rep(FALSE,nrow(genos2)))
              newLocusLocs <- c(mutLocs[k],newLocusLocs)
            }
          }
          colnames(genos1) <- colnames(genoMat1)
          colnames(genos2) <- colnames(genoMat2)
          genos1 <- rbind(genos1,genoMat1)
          genos2 <- rbind(genos2,genoMat2)
        }
        
        newLocusInfo <- polyInfo
        genosIDs <- c(genosIDs,theseOffIDs[kl])
      }
      print(j)
    }
    lethals <- which(newLocusInfo[,4] == "1")

    # prune all individuals that are not currently in the population
    keepIndividualIndices <- which(genosIDs %in% initPopCurrent[,1])
    genos1 <- genos1[keepIndividualIndices,]
    genos2 <- genos2[keepIndividualIndices,]
    genosIDs <- genosIDs[keepIndividualIndices]
    
    # remove loci fixed for the ancestral allele (these do not affect fitness in the population)
    monos <- colSums(rbind(genos1,genos2)) == 0
    if(sum(monos) > 0){
      genos1 <- genos1[,-which(monos == TRUE)]
      genos2 <- genos2[,-which(monos == TRUE)]
      newLocusInfo <- newLocusInfo[-which(monos == TRUE),]
    }
    #--------------------------------------------------------------------------------------------------------
    # prune the loci that are present so that the inbreeding load for survival to 40 years is B ~ 7.26
    #--------------------------------------------------------------------------------------------------------
    sampB <- NULL
    p <- colSums(rbind(genos1,genos2))/(2*nrow(genos1))
    if(fullRecessive == FALSE) hCoef <- 0.5*(exp(-beta*as.numeric(newLocusInfo[,4])))
    if(fullRecessive == TRUE) hCoef <- rep(0,nrow(newLocusInfo))
    singleLocusB <- p*as.numeric(newLocusInfo[,4]) - (p^2)*as.numeric(newLocusInfo[,4]) - 2*(p*(1-p)*as.numeric(newLocusInfo[,4])*hCoef)
    sampB<- sum(p*as.numeric(newLocusInfo[,4])) - sum((p^2)*as.numeric(newLocusInfo[,4])) - 2*sum(p*(1-p)*as.numeric(newLocusInfo[,4])*hCoef)
    
    
    singleLocusB2 <- singleLocusB
    newLocusInfo2 <- newLocusInfo
    bVec <- NULL
    nLocusVec <- NULL
    if(sampB >= B){
      testThisB <- FALSE
      while(testThisB == FALSE){
        if(sum(singleLocusB2) >= B){
          candLocus <- sample(1:nrow(newLocusInfo2),1,replace=FALSE)
          newLocusInfo2 <- newLocusInfo2[-candLocus,]
          singleLocusB2 <- singleLocusB2[-candLocus]
        }
        if(sum(singleLocusB2) <= B){
          testThisB <- TRUE
        }
        bVec <- c(bVec,sum(singleLocusB2))
        nLocusVec <- c(nLocusVec,length(singleLocusB2))
      }
      keepLoci <- match(newLocusInfo2[,2],newLocusInfo[,2])
      
      genos1 <- genos1[,keepLoci]
      genos2 <- genos2[,keepLoci]
      newLocusInfo <- newLocusInfo[keepLoci,]
    }
    
    if(fullRecessive == FALSE) hCoef <- 0.5*(exp(-beta*as.numeric(newLocusInfo[,4])))
    if(fullRecessive == TRUE) hCoef <- rep(0,nrow(newLocusInfo))
    
    #---------------------------------------
    # calculate the current inbreeding load
    #---------------------------------------
    BVec <- NULL
    p <- colSums(rbind(genos1,genos2))/(2*nrow(genos1))
    if(fullRecessive == FALSE) hCoef <- 0.5*(exp(-beta*as.numeric(newLocusInfo[,4])))
    if(fullRecessive == TRUE) hCoef <- rep(0,nrow(newLocusInfo))
    BVec <- c(BVec,sum(p*as.numeric(newLocusInfo[,4])) - sum((p^2)*as.numeric(newLocusInfo[,4])) - 2*sum(p*(1-p)*as.numeric(newLocusInfo[,4])*hCoef))
    hets <- rowSums((genos1 != genos2))/ncol(genos1)
 
    # calculate survival probability
    fitMat <- matrix(NA,nrow=nrow(genos1),ncol= ncol(genos1))   # store the locus specific fitness values
    fitVec <- rep(1,nrow(genos1))                               # save raw relative individual fitness (survival to 40)
    fitVec2 <- rep(1,nrow(genos1))                               # save raw relative individual fitness (annual survival)
    for(j in 1:nrow(genos1)){
      altCounts <- (as.numeric(genos1[j,]) + as.numeric(genos2[j,]))       # count alternative alleles at each locus
      fitPenaltyVec <- rep(0,length(altCounts))                            # locus-specific fitness reduction
      fitPenaltyVec[which(altCounts == 1)] <- (as.numeric(newLocusInfo[which(altCounts == 1),4])) * hCoef[which(altCounts == 1)] # fitness reduction for heterozygous loci
      fitPenaltyVec[which(altCounts == 2)] <- as.numeric(newLocusInfo[which(altCounts == 2),4])                                  # fitness reduction for homozygous loci 
      fitMat[j,] <- fitPenaltyVec
      if(sum(fitPenaltyVec > 0) > 0){
        delLoci <- fitPenaltyVec[which(fitPenaltyVec > 0)]         # remove loci with no derived alleles as the ancestral allele does not reduce fitness
        fitVec[j] <- prod(1-delLoci)                               # calculate the multiplicative fitness reduction due to deleterious genotypes
        fitVec2[j] <- prod(1-delLoci)^(1/39)   
      }
    }
    
    #-----------------------------------------------------------------------------
    # get survival to 40 probability using the regression intercept as reference
    #-----------------------------------------------------------------
    # calculate reference survival to 40 for individual that made it through the first year of life, to be consistent with our sampling of only individuals that made it through very early life
    survProbRaw <- prod(rep(0.9934952,39))    # 0.9956314 is the yearly survival probability of a minimally inbred adult female. Note that prod(rep(0.9956314,39)) == 0.8430335, which is nearly identical to the inferred probability of a non-inbred female surviving to 40 (see supplementary demographic figures)
    Fh <- (max(hets)-hets)/max(hets)                # convert heterozygosity to inbreeding under the assumption that the most heterozygous individual is non-inbred: H = H0(1-F), see Kardos et al. 2016, Evolutionary Applications
    mod <- lm(fitVec~Fh)
    refFit <- mod[[1]][1]                  # intercept of a linear model of fitness ~ inbreeding is the reference fitness
    adjFit <- (fitVec/refFit)*survProbRaw           # adjust raw fitnesses so that the most heterozygous individual has the survival probability associated with the most heterozygous female in our empirical SRKW data
    plot(Fh,adjFit)
    outProbMat <- as.data.frame(cbind(adjFit,Fh))                   
    plot(outProbMat$Fh,outProbMat$adjFit,xlab="F",ylab="survival probability")
    lines(lowess(outProbMat$Fh,outProbMat$adjFit,f=1/3))
    
    #----------------------------------
    # get annual survival probability
    #----------------------------------
    # calculate reference annual survival for individual that made it through the first year of life, to be consistent with our sampling of only individuals that made it through very early life
    survProbRaw <- 0.9934952    # 0.9934952 is the yearly survival probability of a minimally inbred adult male. Note that prod(rep(0.9934952,39)) = 0.7752917, which is nearly identical to the inferred probability of a non-inbred male surviving to 40 (see supplementary demographic figures)
    Fh <- (max(hets)-hets)/max(hets)                # convert heterozygosity to inbreeding under the assumption that the most heterozygous individual is non-inbred: H = H0(1-F), see Kardos et al. 2016, Evolutionary Applications
    mod <- lm(fitVec2~Fh)
    plot(Fh,fitVec2)
    refFit2 <- mod[[1]][1]                  # intercept of a linear model of fitness ~ inbreeding is the reference fitness
    adjFit2 <- (fitVec2/refFit2)*survProbRaw           # adjust raw fitness so that the most heterozygous individual has the survival probability associated with the most heterozygous female in our empirical SRKW data
    plot(Fh,adjFit2)
    outProbMat <- as.data.frame(cbind(adjFit2,Fh))                   
    plot(outProbMat$Fh,outProbMat$adjFit,xlab="F",ylab="survival probability")
    lines(lowess(outProbMat$Fh,outProbMat$adjFit2,f=1/5))
    
    #----------------------------------------------------------------------------------
    # calculate individual genetic loads using maximum intrinsic fitness as reference
    #----------------------------------------------------------------------------------
    if(fullRecessive == FALSE)getIndivLoads(genomeInfo=newLocusInfo,inGenomes1=genos1,inGenomes2=genos2,expRate=beta,fullRecessive=FALSE,multModel=TRUE)
    if(fullRecessive == TRUE)getIndivLoads(genomeInfo=newLocusInfo,inGenomes1=genos1,inGenomes2=genos2,expRate=beta,fullRecessive=TRUE,multModel=TRUE)
    hets <- rowSums((genos1 != genos2))/ncol(genos1)
    fitVec3 <- 1-indivLoads                               # One minus individual genetic load is the metric of expected relative fitness
    Fh <- (max(hets)-hets)/max(hets)                      # convert heterozygosity to inbreeding under the assumption that the most heterozygous individual is non-inbred: H = H0(1-F), see Kardos et al. 2016, Evolutionary Applications
    
    #----------------------------------------------------
    # calculate reference fitness from the genomic data
    #----------------------------------------------------
    refFit3 <- max(fitVec3)                  # max relative fitness in the reference year one is the reference
  
 
    #---------------------------------------------------
    # calculate annual survival under inbreeding depression
    #---------------------------------------------------
    adjFit3 <-(fitVec3/refFit3) * survProbRaw        # fitVec2/refFit2 is expected relative fitness. Raising (fitVec2/refFit2) to (1/39) converts the genetic effects on 
    # fitness to act on yearly survival (we initially simulate genetic load on survival to 40, which is closer to total 
    # fitness than yearly survival). We then multiply the result times baseline age- and sex-specific annual 
    # survival to account for inbreeding effects

    ########
    ######## Save output
    ########
    realizedB [y] <- -1*(lm(log(adjFit)~Fh)[[1]][2])            # calculate lethal equivalents using the regression method
    outSurvProbs[y,] <- adjFit
    outSurvProbs2[y,] <- adjFit2
    outFh[y,] <- Fh
    outHet [y,] <- hets
    print(y)
}   



plot(c(0,0.6),c(0,1.5),xlab="F",ylab="survival probability",type="n")
library(scales)
for(i in 1:nrow(outSurvProbs)){
  lines(lowess(outFh[i,],outSurvProbs[i,]),col=alpha("blue",alpha=0.4))
}

write.table(outSurvProbs,file=paste("survTo40Probs_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,"_sex_",SimSex,sep=""),quote=FALSE,row.names=FALSE)
write.table(outSurvProbs2,file=paste("annSurvProbs_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,"_sex_",SimSex,sep=""),quote=FALSE,row.names=FALSE)
write.table(outFh,file=paste("Fh_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,"_sex_",SimSex,sep=""),quote=FALSE,row.names=FALSE)
write.table(outHet,file=paste("het_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,"_sex_",SimSex,sep=""),quote=FALSE,row.names=FALSE)

