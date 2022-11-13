# Numerically estimate the number of generations needed for a hypothetically completely isolated
# SRKW population to lose the amount of genetic variation equal to the difference
# in heterozygosity between the ARKW and SRKW
setwd("/Users/martin.kardos/Documents/orca/analyses_31May2022/Ne/GONE/lowMissingIndivs/srkw")
srkwHet <- 0.0002877839
arkwHet <- 0.00034646
hetDiff <- arkwHet - srkwHet
NeVec <- read.table("SRKW_pastNe",header=TRUE)[,1]

H0 <- arkwHet
Ht <- srkwHet

# assume GONE-estimated historical Ne
hetIter <- H0
i <- 0
while(hetIter > Ht){
  hetIter <- hetIter - (1/(2*NeVec[i+1]))*hetIter
  print(hetIter)
  print(i)
  i <- i + 1
}
gens_GONE <- i-1
print(gens_GONE)



# assume historical Ne was equivalent to contemporary Ne
hetIter <- H0
i <- 0
while(hetIter > Ht){
  hetIter <- hetIter - (1/(2*27.4))*hetIter
  print(hetIter)
  print(i)
  i <- i + 1
}
gens_constNe <- i-1
print(gens_constNe)
