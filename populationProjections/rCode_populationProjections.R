# cleanup workspace
rm(list=ls())
gc()

setwd("~/Documents/orca/analyses_31May2022/Projections")
# survival model is truncated for younger animals.
# input mutation parameters
Shape <- 0.2
Scale=0.1
PropLethal <- 0.01          # proportion of lethal mutations
beta <- 50                  # exponential rate of decreasing dominance coefficient (h) with increasing size of the selection coefficient (s). 
                            # The dominance model is from Deng & Lynch (1996, Genetics 144 349-360).
B <-  3.24                  # haploid lethal equivalents for survival to 40 years (the empirical sex-averaged value for survival to 40 years)
reps <- 500                 # how manys simulation replicates?
simLethalAlleles <- TRUE    # set to TRUE to leave in lethal alleles, set to FALSE to remove them
NYEARS =  200               # number of years to project
KVec <- rep(200,NYEARS)     # vector of carrying capacities. We used this to avoid runaway population growth in the simulations without inbreeding depression
fullRecessive <- FALSE      # run the simulations using fully recessive deleterious alleles?
simInbDep  <- TRUE          # do you want to include inbreeding depression? Say no if you want to see how the population would grow if everyone had survival as high as the least inbred SRKW.

library(quantPop)
library(pracma)
#--------------------------------------------------------------------------
# read the input genomes from the source population
#--------------------------------------------------------------------------
inGenos1 <- read.table(paste("inputKillerGenos1_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,sep=""),header=FALSE)
inGenos2 <- read.table(paste("inputKillerGenos2_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,sep=""),header=FALSE)
inLocusInfo <- read.table(paste("inputKillerLocusInfo_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,sep=""),header=FALSE)
colnames(inGenos1) <- inLocusInfo[,2]
colnames(inGenos2) <- inLocusInfo[,2]

# remove lethals and set the lethal mutation rate to zero if necessary
if(simLethalAlleles == FALSE){
  PropLethal <- 0.0    # switch to proportion of lethals = 0 for the simulations below
  lethals <- which(inLocusInfo[,4] == "1")
  inGenos1 <- inGenos1[,-lethals]
  inGenos2 <- inGenos2[,-lethals]
  inLocusInfo <- inLocusInfo[-lethals,]
}

#-------------------------------------------------------------------------------------------------------------------
# remove any individuals that are homozygous for a recessive lethal. By definition they cannot be pedigree founders
# because they would be dead
#-------------------------------------------------------------------------------------------------------------------
# remove any individuals that are homozygous for lethal recessives
#identify loci with lethal alleles
if(sum(inLocusInfo[,4] == 1) > 0){      # if lethal alleles are present
  lethals <- which(inLocusInfo[,4] ==1)
  # see if any lethals are homozygous in any individual
  homLethals <- colSums((inGenos1[,lethals] == TRUE) & (inGenos2[,lethals] == TRUE))
  if(sum(homLethals) > 0){
    lethalHomLoci <- lethals[which(homLethals > 0)]
    homLethalInds <- NULL
    for(i in 1:length(lethalHomLoci)){
      homLethalInds <- c(homLethalInds,which(inGenos1[,lethalHomLoci[i]] == TRUE & inGenos2[,lethalHomLoci[i]]))
    }
    homLethalInds <- sort(unique(homLethalInds))
   inGenos1 <- inGenos1[-homLethalInds,]
   inGenos2 <- inGenos2[-homLethalInds,]
  }
}
  
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

survivModel <- readRDS("fit.rds")
# load in demographic data and parameters
load("demog_models_projection.Rdata")
years_to_sample = 2000                       # assume the envirnomental component is fixed and equivalent to the that of the year 2000
pFemale = 50/100                             # sex ratio at birth is even
sdFemale = sqrt(pFemale*(1-pFemale)/100)     # stochastic component of sex ratio
currentPop$froh = rep(-1.2622,nrow(currentPop)) # set froh to the minimum standardized (substract mean and divide by standard deviation) value observed

#store outputs
loadMat<- matrix(NA,nrow=reps,ncol=NYEARS+1)               # total genetic load
meanFMat <- matrix(NA,nrow=reps,ncol=NYEARS+1)             # mean pedigree-based inbreeding through time
meanHMat <- matrix(NA,nrow=reps,ncol=NYEARS+1)             # mean heterozygosity through time
popSize = array(0, dim=c(reps, NYEARS))                    # population size through time
BMat <- matrix(NA,nrow=reps,ncol=NYEARS+1)                 # inbreeding load through time

# initPopCurrent keeps track of individuals present in the population through time
initPopCurrent = currentPop[,c("animal", "pod", "matriline", "sexF1M2", "age","froh")] # animal, pod, matriline, sex, age
initPopCurrent = dplyr::rename(initPopCurrent, sex=sexF1M2)                # coding of sex
initPopCurrent$animal = as.character(initPopCurrent$animal)                # individual ID
initPopCurrent$matriline = as.numeric(as.factor(initPopCurrent$matriline)) # which matriline an individual is in
initPopCurrent$dad = NA     # keep track of dads



numWhale = dim(initPopCurrent)[1]
initPopCurrent$hasCalf = 0
initPopCurrent$hasCalf[which(as.character(currentPop$animal) %in% as.character(currentPop$mom[which(currentPop$age==1)]))] = 1
initPopCurrent$mom = as.character(currentPop$mom)

# assign sexes to unsexed individuals [0s]
initPopCurrent$sex = as.numeric(initPopCurrent$sex)
initPopCurrent$sex[which(initPopCurrent$sex==0)] = ifelse(runif(length(which(initPopCurrent$sex==0))) < rnorm(length(which(initPopCurrent$sex==0)),pFemale, sdFemale), 1, 2)
initPopCurrent$year <- 1
initPopCurrent$sexF1M2 <- initPopCurrent$sex
newID = 9999         # start iterating new whale IDs

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

# update the pedigree
newPed <- matrix(NA,nrow=length(onlyMomIDs) + length(missingMoms),ncol=ncol(SRPed))
newPed[,1] <- c(onlyMomIDs,missingMoms)
newPed[1:length(onlyMomIDs),6] <- initPopCurrent[match(onlyMomIDs,initPopCurrent[,1]),9]   # moms of individuals with only mom ids
newPed[1:length(onlyMomIDs),2] <- c(rep("2020",4),"1985",rep("2020",6))                    # say that they are born in 2020 to ensure the moms are alive when simulating their genomes
newPed[1:length(onlyMomIDs),4] <- initPopCurrent[match(onlyMomIDs,initPopCurrent[,1]),4]   # sex of individuals with only mom IDs available
newPed[(length(onlyMomIDs)+1):(nrow(newPed)),4] <- 1                                       # moms are all female
newPed <- as.data.frame(newPed)
colnames(newPed) <- colnames(SRPed)
SRPed <- rbind(SRPed,newPed)

# add founders to accomodate mating to produce individuals with only one known parent. The other parent is randomly drawn from the source population.
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

#store outputs
meanFMat <- matrix(NA,nrow=reps,ncol=NYEARS+1)
meanHMat <- matrix(NA,nrow=reps,ncol=NYEARS+1)
popSize = array(0, dim=c(reps, NYEARS))
BMat <- matrix(NA,nrow=reps,ncol=NYEARS+1)
for(y in 1:reps){
  initPopCurrent2 <- initPopCurrent
  #---------------------------------------------------------------------------------------------------
  # randomly assign each pedigree founder a genome from the source population
  #---------------------------------------------------------------------------------------------------
  founders <- SRPed[is.na(SRPed$dad) & is.na(SRPed$mom),]
  nonFounders <- SRPed[which((is.na(SRPed$dad) & is.na(SRPed$mom)) == FALSE),]
  sampGens <- sort(sample(1:nrow(inGenos1),nrow(founders),replace=FALSE))    # which genomes to assign to the pedigree founders
  newLocusInfo <- as.matrix(inLocusInfo)
  newLocusInfo[,1] <- rep(0,nrow(newLocusInfo))
  genos1 <- inGenos1[sampGens,]                               # first haploid genome copies
  genos2 <- inGenos2[sampGens,]                               # second haploid genome copies
  colnames(genos1) <- as.character(newLocusInfo[,2])
  colnames(genos2) <- as.character(newLocusInfo[,2])
  
  # add founder genomes that will mate with a known parent to produce offspring with only one known parent
  genosIDs <- founders[,1]                              # link IDs to genomes

  #--------------------------------------------------------------------------------------------------------------------------
  # prune loci that are present in the founders so that the initial inbreeding load for survival to 40 years is B
  # This is to make simulating the pedigree faster/
  #--------------------------------------------------------------------------------------------------------------------------
  B0 <- B
  sampB <- NULL
  p <- colSums(rbind(genos1,genos2))/(2*nrow(genos1))                                      # allele frequencies
  if(fullRecessive == FALSE) hCoef <- 0.5*(exp(-beta*as.numeric(newLocusInfo[,4])))        # set the dominance coefficients according to beta
  if(fullRecessive == TRUE) hCoef <- rep(0,nrow(newLocusInfo))
  singleLocusB <- p*as.numeric(newLocusInfo[,4]) - (p^2)*as.numeric(newLocusInfo[,4]) - 2*(p*(1-p)*as.numeric(newLocusInfo[,4])*hCoef)       # single locus contributions to inbreeding load
  sampB<- sum(p*as.numeric(newLocusInfo[,4])) - sum((p^2)*as.numeric(newLocusInfo[,4])) - 2*sum(p*(1-p)*as.numeric(newLocusInfo[,4])*hCoef)  # total inbreeding load
  
  
  singleLocusB2 <- singleLocusB
  newLocusInfo2 <- newLocusInfo
  bVec <- NULL
  nLocusVec <- NULL
  if(sampB >= B0){
    testThisB <- FALSE
    while(testThisB == FALSE){
      if(sum(singleLocusB2) >= B0){                                     # randomly identify a locus to remove
        candLocus <- sample(1:nrow(newLocusInfo2),1,replace=FALSE)
        newLocusInfo2 <- newLocusInfo2[-candLocus,]
        singleLocusB2 <- singleLocusB2[-candLocus]
      }
      if(sum(singleLocusB2) <= B0){                                     # see if the current inbreeding load is lower than B
        testThisB <- TRUE
      }
      bVec <- c(bVec,sum(singleLocusB2))
      nLocusVec <- c(nLocusVec,length(singleLocusB2))
      #print(paste("pruned load = ",sum(singleLocusB2) ))              # uncomment if you want to see the progression of the pruning process
    }
    keepLoci <- match(newLocusInfo2[,2],newLocusInfo[,2])              # final set of loci to start with  
    genos1 <- genos1[,keepLoci]
    genos2 <- genos2[,keepLoci]
    newLocusInfo <- newLocusInfo[keepLoci,]
  }
  hets <- rowSums((genos1 != genos2))/ncol(genos1)                     # individual heterozygosity
  if(fullRecessive == FALSE) getIndivLoads(genomeInfo=newLocusInfo,inGenomes1=genos1,inGenomes2=genos2,expRate=beta,fullRecessive=FALSE,multModel=TRUE)
  if(fullRecessive == TRUE) getIndivLoads(genomeInfo=newLocusInfo,inGenomes1=genos1,inGenomes2=genos2,expRate=beta,fullRecessive=TRUE,multModel=TRUE)
  
  # set the reference genetic load as that expected for individual with the highest heterozygosity (assumed to be nob-inbred). This works because we are modelling inbreeding depression within the population. Drift load is assumed fixed and reflected in the fitness of all individuals in the population.
  Fh <- (max(hets)-hets)/max(hets)        # inbreeding defined as the proportional reduction in heterozygosity compared to the most heterozygous (least inbred) individual
  refFit <- max(1-indivLoads)             # the reference fitness
  
  #---------------------------------------------------
  # simulated mating in the known pedigree
  #---------------------------------------------------
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
        # make a candidate offspring genome. It has to survive to breeding age to be included in the pedigree
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
        fitVec <- 1-indivLoads            # intrinsic fitness
        newDF <- NULL                     # data frame to store demographic information
        newDF$sexF1M2 <- rep(as.numeric(theseSexes[kl]),10)
        newDF$age <- 1:10
        newDF$froh <- rep(-1.2622,10)
        newDF$year <- rep(2000,10)
        p_surv = fitted(surv_mod,newDF)
        survProbRaw <- 0.5*prod(p_surv[,1])                  # assume baseline survival of neonates is 0.5
        adjSurvProb <-((fitVec/refFit))^(1/39) * survProbRaw        # fitVec2/refFit2 is expected relative fitness. Raising (fitVec2/refFit2) to (1/39) converts the genetic effects on 
        
        # incorporate inbreeding depression on survival to reproductive age
        survTest <- runif(1,min=0,max=1) <= adjSurvProb      # does the simulated individual survive to 10 years?
        if(survTest == TRUE)survived <- 1                    # continue if the offspring survived
        if(simInbDep == FALSE) survived <- 1                 # make the offspring survive no matter the genome if you're not simulating inbreeding depression
      }
      
      if(ncol(genoMat1) == ncol(genos1)){
        genos1 <- rbind(genos1,genoMat1)          # update genotype matrices and locus information with new individuals
        genos2 <- rbind(genos2,genoMat2)
      }
      
      # add new mutations if there were any in the simulated offspring
      if(ncol(genoMat1) != ncol(genos1)){                        
        mutLocs <- as.numeric(polyInfo[which(polyInfo[,6] %in% newLocusInfo[,6] == FALSE),6])        # mutation locations
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
  keepIndividualIndices <- which(genosIDs %in% initPopCurrent2[,1])
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
  
  #---------------------------------------------------------------------------------------------------------------------------------------------------------------
  # prune the loci that are present so that the inbreeding load for survival to 40 years is B, in case it went up during the simulation of the current pedigree
  #---------------------------------------------------------------------------------------------------------------------------------------------------------------
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
      if(sum(singleLocusB2) >= B){                                     # randomly selection a locus to remove
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
  
  #----------------------------------------
  # initialize pedigree information
  #----------------------------------------
  foundPed <- SRPed[which(is.na(SRPed[,6])),c(1,6,5)]
  nonFoundPed <- SRPed[which(is.na(SRPed[,6]) == FALSE),c(1,6,5)]
  outPed <- rbind(foundPed,nonFoundPed)
  outPed[,1] <- as.character(outPed[,1])
  ped <- as.matrix(outPed)
  colnames(ped) <- c("ID","mom","dad")
  pedList <- list() # store an updated pedigree each year
  pedList[[1]] <- ped

  #iterate years and initialize checks for extinction
  # initialize yrs
  sizeNow <- dim(initPopCurrent2)[1]   # current population size
  bothSexes <- TRUE                   # are both sexes present?
  meanHVec <- NULL                    # store heterozygosity from this repetition
  meanLoadVec <- NULL                 # store individual genetic loads 
  BVec <- NULL                        # store inbreeding load
  yrs <- 1
  while(yrs<= NYEARS & sizeNow > 2 & bothSexes == TRUE) {        # stop if/when population size hits <=2 or only one sex is left
    print(paste("year =",yrs))
    
    # first step is find females available to give birth
    indx = which(initPopCurrent2$sex == 1 & as.numeric(initPopCurrent2$age) >= 10 & as.numeric(initPopCurrent2$age) < 43 & initPopCurrent2$hasCalf == 0)
    
    # second step is to calculate predicted fecundity rates and make fecundity stochastic - every individual's pregnancy is a coin flip
    # if a range of years is used to draw demographic rates from sample those randomly
    initPopCurrent2$year = as.numeric(sample(c(as.character(years_to_sample)), size=1))
    if(length(indx) > 0) {
      # bind together the current moms with the matching fecundity data
      p_fec = predict(fec_mod, initPopCurrent2[indx,])[,1]
      pregMoms <- NULL
      realized_fec <- runif(length(p_fec)) < p_fec
      if(sum(realized_fec) > 0) pregMoms <- indx[which(realized_fec == TRUE)]

      # third step is to make moms that aren't mate limited give birth to calves of known sexes
      theseMoms <- NULL
      theseDads <- NULL
      newOffIDs <- NULL
      if(length(pregMoms) > 0) {
        thisPed <- NULL # pedigree information for this year
        for(ll in 1:length(pregMoms)) {
          # loop over moms and only let them breed if there's a mature male in a different matriline
          dads = which(initPopCurrent2$sex == 2 & as.numeric(initPopCurrent2$age) > 12)
          if(length(dads) > 0) {
            
            # assign the pod / matriline to be the same of the mom
            newpod = initPopCurrent2$pod[pregMoms[ll]]
            newmat = initPopCurrent2$matriline[pregMoms[ll]]
            
            # sex is stochastic
            newsex = ifelse(runif(1) < rnorm(1, pFemale, sdFemale), 1, 2)
            
            # bookkeeping
            newage = 0
            newcalf = 0
            newmom = initPopCurrent2$animal[pregMoms[ll]]
            theseMoms <- c(theseMoms,newmom)
            
            # sample from potential dads in proportion to their estimated relative reproductive output
            # their ages are initPopCurrent2$age[dads], and the fecundity model defined above outside of loops
            inflation_factor = 5
            pred.male.fecundity[32:200] = pred.male.fecundity[31]    # set fecundity expectation for old males
            probs = pred.male.fecundity[as.numeric(initPopCurrent2$age[dads])]
            probs[which.max(initPopCurrent2$age[dads])] = probs[which.max(initPopCurrent2$age[dads])] * inflation_factor
            if(length(dads) > 1) newdad = initPopCurrent2$animal[sample(dads,1, prob=probs)] # sample dads and allow for only one potential father
            if(length(dads) == 1) newdad = initPopCurrent2$animal[dads]
            theseDads <- c(theseDads,newdad)
            
            # iterate offspring IDs
            newID = newID + 1
            newOffIDs <- c(newOffIDs,newID)
            
            # add calves to the population
            newdf = data.frame(animal=newID,
                               pod=newpod,matriline=newmat,sex=newsex,
                               age=newage,froh=-1.2622,dad=newdad,hasCalf=newcalf,
                               mom=newmom,year=initPopCurrent2$year[1],sexF1M2=theseSexes[ll])
            if(colnames(initPopCurrent2)[ncol(initPopCurrent2)] == "year") initPopCurrent2$sexF1M2 <- initPopCurrent2$sex
            initPopCurrent2 = rbind(initPopCurrent2, newdf)
            thisPed <- rbind(thisPed,c(newID,theseMoms[ll],theseDads[ll]))
          }
        }
        ped <- rbind(ped,thisPed)      # update the pedigree
        
        #----------------------------------------
        # meiosis and mutation
        #----------------------------------------
        mutRate <- 5e-8                                    # one deleterious mutation per diploid offspring
        if(is.null(newOffIDs) == FALSE){
          mating_delMutation(delMu <- mutRate,
                             genSize=10000000,
                             chroms=20,
                             mothers= theseMoms,
                             fathers= theseDads,
                             genoIDVec = genosIDs,
                             propLethal=PropLethal,
                             delMuGammaShape=Shape,
                             delMuGammaScale=Scale,
                             Beta=beta,
                             importGenos=TRUE,
                             genoMat1Name=genos1,
                             genoMat2Name=genos2,
                             inLocusInfo=newLocusInfo,
                             year=yrs)
          
          #-------------------------------------------------------------------
          # update the population genomes with new offspring and new mutations
          #-------------------------------------------------------------------
          if(ncol(genoMat1) > ncol(genos1)){    # add loci AND individuals if there were any mutations in the new offspring
            newLociLocs <- polyInfo[which(polyInfo[,2] %in% newLocusInfo[,2] == FALSE),6]
            oldLociLocs <- newLocusInfo[,6]
            combineLocs <- as.numeric(c(oldLociLocs,newLociLocs))
            
            newGenos1 <- genos1   # new genotype matrices including newly mutated loci in the new offspring
            newGenos2 <- genos2
            newGenos1 <- cbind(newGenos1,matrix(FALSE,nrow=nrow(newGenos1),ncol=length(newLociLocs)))
            newGenos2 <- cbind(newGenos2,matrix(FALSE,nrow=nrow(newGenos2),ncol=length(newLociLocs)))
            
            # sort loci by genomic location to conform with the new genomes
            newGenos1 <- newGenos1[,order(combineLocs)]
            newGenos2 <- newGenos2[,order(combineLocs)]
            
            colnames(newGenos1) <- colnames(genoMat1)
            colnames(newGenos2) <- colnames(genoMat2)
            # add the new genomes to the old
            newGenos1 <- rbind(newGenos1,genoMat1)
            newGenos2 <- rbind(newGenos2,genoMat2)
          }
          if(ncol(genoMat1) == ncol(genos1)){    # add loci AND individuals if there were any mutations in the new offspring
            # new genotype matrices
            newGenos1 <- genos1
            newGenos2 <- genos2
            colnames(newGenos1) <- colnames(genoMat1)
            colnames(newGenos2) <- colnames(genoMat2)
            # add the new genomes to the old
            newGenos1 <- rbind(newGenos1,genoMat1)
            newGenos2 <- rbind(newGenos2,genoMat2)
          }
          # replace the old genomic data with the new
          genos1 <- newGenos1
          genos2 <- newGenos2
          colnames(genos1) <- colnames(genoMat1)
          colnames(genos2) <- colnames(genoMat2)
          newLocusInfo <- polyInfo
          genosIDs <- c(genosIDs,newOffIDs)
        }
      } # end if(length(pregMoms)
      
      if(sum(theseMoms == theseDads) > 0) break
    }
    # update the pedigree time series
    pedList[[yrs+1]] <- ped
    
    # bookkeeping: update whales that have calves
    if(sum(initPopCurrent2$hasCalf==1) > 0) initPopCurrent2$hasCalf[which(initPopCurrent2$hasCalf==1)] = 0
    if(length(pregMoms) > 0 & length(indx) > 0) initPopCurrent2$hasCalf[pregMoms] = 1                  # Marty modified this to also require indx have length > 0
    
    initPopCurrent2$sexF1M2 = as.numeric(initPopCurrent2$sex)
    initPopCurrent2 = dplyr::select(initPopCurrent2, -sexF1M2)
    
    #--------------------------------------------------
    # calculate individual genetic loads
    #--------------------------------------------------
    if(fullRecessive == FALSE)getIndivLoads(genomeInfo=newLocusInfo,inGenomes1=genos1,inGenomes2=genos2,expRate=beta,fullRecessive=FALSE,multModel=TRUE)
    if(fullRecessive == TRUE)getIndivLoads(genomeInfo=newLocusInfo,inGenomes1=genos1,inGenomes2=genos2,expRate=beta,fullRecessive=TRUE,multModel=TRUE)
    hets <- rowSums((genos1 != genos2))/ncol(genos1)
    fitVec2 <- 1-indivLoads                               # One minus individual genetic load is the metric of intrinsic fitness
    Fh <- (max(hets)-hets)/max(hets)                      # convert heterozygosity to inbreeding under the assumption that the most heterozygous individual is non-inbred: H = H0(1-F), see Kardos et al. 2016, Evolutionary Applications
    
    #----------------------------------------------------
    # calculate reference fitness from the genomic data
    #----------------------------------------------------
    if(yrs == 1) refFit2 <- max(fitVec2)                  # max relative fitness in the reference year one is the reference
    
    #--------------------------------------------------------------
    # assign sex- and age-specific baseline survival probability
    #--------------------------------------------------------------
    initPopCurrent2$sexF1M2 <- as.numeric(initPopCurrent2$sex)
    survProbRaw <- fitted(surv_mod,initPopCurrent2)[,1]         # the baseline fitness for each individual is fitted value for sex- and age-specific survival probability associated with minimally inbred individuals in the empirical data
    
    # make baseline survival of all age zeros 0.5, consistent with with empirical data from West Coast N America populations; see Oleksiuk et al. 1990, Rep. Int. Whal. Commn, Issue 12)
    if(sum(initPopCurrent2$age == 0) > 0){
      survProbRaw[which(initPopCurrent2$age == 0)] <- 0.5
    }
    
    #---------------------------------------------------
    # calculate fitness under inbreeding depression
    #---------------------------------------------------
    adjFit2 <-((fitVec2/refFit2)^(1/39)) * survProbRaw        # fitVec2/refFit2 is expected relative fitness. Raising (fitVec2/refFit2) to (1/39) converts the genetic effects on 
                                                              # fitness to act on yearly survival (we initially simulate genetic load on survival to 40, which is closer to total 
                                                              # fitness than yearly survival). We then multiply the result times baseline age- and sex-specific annual 
                                                              # survival to account for inbreeding effects

    #-------------------------------
    # save genetic load information
    #-------------------------------
    meanLoadVec <- c(meanLoadVec,mean(indivLoads))                                         # save individual loads
    if(fullRecessive == FALSE) hCoef <- 0.5*(exp(-beta*as.numeric(newLocusInfo[,4])))      # calculate the dominance coefficients
    if(fullRecessive == TRUE) hCoef <- rep(0,nrow(newLocusInfo))
    p <- colSums(rbind(genos1,genos2))/(2*nrow(genos1))                                    # allele frequency
    BVec <- c(BVec,sum(p*as.numeric(newLocusInfo[,4])) - sum((p^2)*as.numeric(newLocusInfo[,4])) - 2*sum(p*(1-p)*as.numeric(newLocusInfo[,4])*hCoef))
    meanHVec <- c(meanHVec,mean(rowSums(genos1[genosIDs %in% initPopCurrent2[,1],] != genos2[genosIDs %in% initPopCurrent2[,1],])/ncol(genos1)))
    
    #-------------------------------------
    # incorporate inbreeding depression
    #-------------------------------------
    if(simInbDep == TRUE) probSurv <- adjFit2                       # survival probability under inbreeding depression
    if(simInbDep == FALSE) probSurv <- survProbRaw                  # assume baseline survival if you're not simulating inbreeding depression
    
    # penalize survival probability if the population is above the current carrying capacity to keep the population sizes reasonable in the simulation. This is a ceiling model of density dependence.
    if(nrow(initPopCurrent2) > KVec[yrs]){
      expectedSurvivors <- sum(adjFit2)
      neededSurvivalRate <- KVec[yrs]/nrow(initPopCurrent2)
      meanSurv <- mean(adjFit2)
      if(meanSurv > neededSurvivalRate){
        survPenalty <- neededSurvivalRate/meanSurv
        probSurv <-  probSurv*survPenalty
      }
    }
    
    # step 5: stochastic survival to kill whales
    liveOrDie = rep(1,length(probSurv))
    dead = which(runif(probSurv) > probSurv)
    liveOrDie[dead] = 0
    
    # step 6: see if any of these dead animals has a calf - if so, kill the calf too
    if(length(which(liveOrDie == 0)) > 0) {
      for(ll in 1:length(dead)) {
        # kill the calf
        if(is.na(initPopCurrent2$hasCalf[dead[ll]]) == FALSE & initPopCurrent2$hasCalf[dead[ll]] == 1) {
          liveOrDie[which(initPopCurrent2$mom == dead[ll])] = 0
        }
      }
    }
    
    # step 7: bookkeeping at the end of the time step
    # first remove dead animals from the population
    if(length(dead) > 0 ) deadWhales = initPopCurrent2[which(liveOrDie==0),]  ## MIKE ADDITION - list of who is dead
    if(length(dead) > 0) initPopCurrent2 = initPopCurrent2[-which(liveOrDie==0),]
    
    # remove dead individuals from the genotype matrices
    remInds <- (genosIDs %in% initPopCurrent2[,1]) == FALSE
    if(sum(remInds) > 0){
      genos1 <- genos1[-which(remInds == TRUE),]
      genos2 <- genos2[-which(remInds == TRUE),]
      genosIDs <- genosIDs[-which(remInds == TRUE)]
    }
    
    print(paste("number of segregating loci = ", ncol(genos1),sep=""))
    print(paste("Lethal equivalents = ", BVec[length(BVec)],sep=""))
    
    # second, age remaining whales
    initPopCurrent2$age = as.numeric(initPopCurrent2$age) + 1
    
    # third, record pop size to later calculate recovery targets
    popSize[y,yrs] = dim(initPopCurrent2)[1]
    print(paste("population size = ",popSize[y,yrs],sep=""))
    
    #iterate years and check for extinction
    yrs <- yrs + 1
    sizeNow <- dim(initPopCurrent2)[1]
    bothSexes <- length(unique(initPopCurrent2$sex)) > 1
    initPopCurrent2 <- initPopCurrent2[,-ncol(initPopCurrent2)]
  } # this is end of yrs loop
  
  # record pedigree inbreeding coefficients after the bottleneck
  meanFVec <- rep(NA,sum((popSize[y,] != 0)) - 1)
  for(k in 1:length(meanFVec)){
    thisPed <- pedList[[k]]
    kinMat <- kinship(id=thisPed[,1],dadid=thisPed[,2],momid=thisPed[,3])
    thisFVec <- rep(NA,nrow(thisPed))
    for(lll in 1:nrow(thisPed)){
      if(is.na(thisPed[lll,2])) thisFVec[k] <- 0
      if(is.na(thisPed[lll,2]) == FALSE) thisFVec [lll]<- kinMat[which(rownames(kinMat) == thisPed[lll,2]),which(colnames(kinMat) == thisPed[lll,3])]
    }
    thisFVec[is.na(thisFVec)] <- 0
    meanFVec[k] <- mean(thisFVec)
  }
    # record genetic outputs
    loadMat[y,1:length(meanLoadVec)] <- meanLoadVec
    meanFMat[y,1:length(meanFVec)] <- meanFVec
    meanHMat[y,1:length(meanHVec)] <- meanHVec
    BMat[y,1:length(BVec)] <- BVec
    print(paste("************************************** simulation ",y," done",sep=""))
}   

if(simInbDep == TRUE){
  write.table(popSize,file=paste("popSize_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,sep=""),row.names=FALSE,quote=FALSE)
  write.table(meanFMat,file=paste("meanF_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,sep=""),row.names=FALSE,quote=FALSE)
  write.table(meanHMat,file=paste("meanH_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,sep=""),row.names=FALSE,quote=FALSE)
  write.table(BMat,file=paste("B_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,sep=""),row.names=FALSE,quote=FALSE)
}

if(simInbDep == FALSE){
  write.table(popSize,file="popSize_noInbDep",row.names=FALSE,quote=FALSE)
  write.table(meanFMat,file="meanF_noInbDep",row.names=FALSE,quote=FALSE)
  write.table(meanHMat,file="meanH_noInbDep",row.names=FALSE,quote=FALSE)
  write.table(BMat,file="B_noInbDep",row.names=FALSE,quote=FALSE)
}
popSize <- read.table(paste("popSize_shape",Shape,"_scale",Scale,"_propLethal",PropLethal,"_beta",beta,sep=""),header=TRUE)
# plot results
library(scales)
newPopSize <- popSize[which(popSize[,10] >0),]
newPopSize <- cbind(rep(74,nrow(newPopSize)),newPopSize)
plot(c(0,200),c(0,90),xlab="Year",ylab="Population Size",type="n")
for(i in 1:nrow(newPopSize)){
  lines(0:200,as.numeric(newPopSize[i,]),col=alpha("darkred",alpha=0.05))
}

means <- colMeans(newPopSize)
quant95  <- matrix(NA,nrow=2,ncol=length(means))
quant50  <- matrix(NA,nrow=2,ncol=length(means))

for(i in 1:ncol(newPopSize)){
  quant95[,i] <- quantile(newPopSize[,i],probs=c(0.025,0.975))
  quant50[,i] <- quantile(newPopSize[,i],probs=c(0.25,0.75))
}
polygon(c(0:200,rev(0:200)),c(quant95[1,],rev(quant95[2,])),col=alpha("darkred",alpha=0.1),border=NA)
polygon(c(0:200,rev(0:200)),c(quant50[1,],rev(quant50[2,])),col=alpha("darkred",alpha=0.3),border=NA)
lines(0:200,means,col="darkred",lwd=3)


