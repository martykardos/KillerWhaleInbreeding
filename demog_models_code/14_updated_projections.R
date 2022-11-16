
# load in demographic data and parameters
load("demog_models_projection.Rdata")
# scenario used affects which rates are used for projections
years_to_sample = 2000

# These were the sex ratios we'd used during the workshops -- this is
# a random variable, rather than fixed, so the uncertainty is propagated
pFemale = 50/100 # this is the number of female births/total up through 2018
sdFemale = sqrt(pFemale*(1-pFemale)/100)

nSims = 1000 # number of simulations
NYEARS = 25 # number of years to project

popSize = array(0, dim=c(nSims, NYEARS))
#currentPop = readRDS("data/current_pop_2021.rds")
currentPop$froh[which(is.na(currentPop$froh))] = rnorm(length(which(is.na(currentPop$froh))),0,1) # in reality, raw values were standardized using mean = 0.03234788, sd = 0.02563815

for(i in 1:nSims) {
  #print(i)
  # get the current age / sex / pod structure of the population
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
  initPopCurrent$sex[which(initPopCurrent$Sex==0)] = ifelse(runif(length(which(initPopCurrent$sex==0))) < rnorm(length(which(initPopCurrent$sex==0)),pFemale, sdFemale), 1, 2)
  newID = 9999
  
  for(yrs in 1:dim(popSize)[2]) {
    # first step is find females available to give birth
    indx = which(initPopCurrent$sex == 1 & as.numeric(initPopCurrent$age) >= 10 & as.numeric(initPopCurrent$age) < 43 & initPopCurrent$hasCalf == 0)
    # second step is to calculate predicted fecundity rates and make fecundity stochastic - every individual's pregnancy is a coin flip
    
    # if a range of years is used to draw demographic rates from sample those randomly
    initPopCurrent$year = as.numeric(sample(c(as.character(years_to_sample)), size=1))
    
    if(length(indx) > 0) {
      # bind together the current moms with the matching fecundity data
      p_fec = predict(fec_mod, initPopCurrent[indx,])[,1]

      pregMoms = indx[which(runif(length(p_fec)) < p_fec)]
      # third step is to make moms that aren't mate limited give birth to calves of known sexes
      if(length(pregMoms) > 0) {
        
        for(ll in 1:length(pregMoms)) {
          # loop over moms and only let them breed if there's a mature male in a different matriline
          dads = which(initPopCurrent$sex == 2 & as.numeric(initPopCurrent$age) > 12)
          #dads = which(initPopCurrent$Sex == 2 & initPopCurrent$breedingGroup != initPopCurrent$breedingGroup[pregMoms[ll]] & as.numeric(initPopCurrent$age1) > 12)
          if(length(dads) > 0) {
            # assign the pod / matriline to be the same of the mom
            newpod = initPopCurrent$pod[pregMoms[ll]]
            newmat = initPopCurrent$matriline[pregMoms[ll]]
            # sex is stochastic
            newsex = ifelse(runif(1) < rnorm(1, pFemale, sdFemale), 1, 2)
            # bookkeeping
            newage = 0
            newcalf = 0
            newmom = initPopCurrent$animal[pregMoms[ll]]
            
            # sample from potential dads in proprtion to their estimated relative reproductive output
            # their ages are initPopCurrent$age1[dads], and the fecundity model defined above outside of loops
            # sample from potential dads in proprtion to their estimated relative reproductive output
            # their ages are initPopCurrent$age1[dads], and the fecundity model defined above outside of loops
            inflation_factor = 5
            pred.male.fecundity[32:200] = pred.male.fecundity[31]
            probs = pred.male.fecundity[as.numeric(initPopCurrent$age[dads])]
            probs[which.max(initPopCurrent$age1[dads])] = probs[which.max(initPopCurrent$age[dads])] * inflation_factor
            newdad = initPopCurrent$animal[sample(dads,1, prob=probs)]

            newID = newID + 1
            # add calves to the population
            newdf = data.frame(animal=newID,
                               pod=newpod,matriline=newmat,sex=newsex,
                               age=newage,dad=newdad,hasCalf=newcalf,
                               froh=rnorm(1,0,1),
                               mom=newmom,year=initPopCurrent$year[1])
            initPopCurrent = rbind(initPopCurrent, newdf)
            
          }# end if(length(dads) > 0) {
        } # end ll loop
      } # end if(length(pregMoms)
      
    }
    # bookkeeping: update whales that have calves
    initPopCurrent$hasCalf[which(initPopCurrent$hasCalf==1)] = 0
    initPopCurrent$hasCalf[pregMoms] = 1
    
    # step 4 is calculate predicted survival at age
    initPopCurrent$sexF1M2 = as.numeric(initPopCurrent$sex)
    p_surv = predict(surv_mod, initPopCurrent)
    p_surv[which(initPopCurrent$age==0)] = 0.5
    initPopCurrent = dplyr::select(initPopCurrent, -sexF1M2)
    
    # step 5: stochastic survival to kill whales
    liveOrDie = rep(1,length(p_surv))
    dead = which(runif(length(p_surv)) > p_surv)
    liveOrDie[dead] = 0
    
    # step 6: see if any of these dead animals has a calf - if so, kill the calf too
    if(length(which(liveOrDie == 0)) > 0) {
      for(ll in 1:length(dead)) {
        # kill the calf
        if(is.na(initPopCurrent$hasCalf[dead[ll]]) == FALSE & initPopCurrent$hasCalf[dead[ll]] == 1) {
          liveOrDie[which(initPopCurrent$mom == dead[ll])] = 0
        }
      }
    }
    
    # step 7: bookkeeping at the end of the time step
    # first remove dead animals from the population
    if(length(dead) > 0 ) deadWhales = initPopCurrent[which(liveOrDie==0),]  ## MIKE ADDITION - list of who is dead
    if(length(dead) > 0) initPopCurrent = initPopCurrent[-which(liveOrDie==0),]
    # second age remaining whales
    initPopCurrent$age = as.numeric(initPopCurrent$age) + 1
    # third record pop size to later calculate recovery targets
    popSize[i,yrs] = dim(initPopCurrent)[1]
    
  } # this is end of yrs loop
  
}