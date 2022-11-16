library(kwdemog)
library(dplyr)
# genetic data
d = read.csv("data/FswithLHandPed.csv") 
d = dplyr::select(d, 
                  id, Froh_min1Mb,
                  Froh_min5Mb,
                  Froh_min10Mb)
# center genetic data -- necessary for prior 0.2748935, 0.05170584
d$Froh_min1Mb = as.numeric(scale((d$Froh_min1Mb)))#as.numeric(scale(log(d$Froh_min1Mb)))
d$Froh_min5Mb = as.numeric(scale((d$Froh_min5Mb)))
d$Froh_min10Mb = as.numeric(scale(d$Froh_min10Mb))
# demography data
data("orca")
expanded = kwdemog::expand(orca)

demog = dplyr::mutate(expanded, pod_char = substr(animal,1,1),
                      num_char = substr(animal,2,length(animal))) %>%
  dplyr::mutate(id = paste0(pod_char,as.numeric(num_char)))

# join in genetic data
demog = dplyr::left_join(demog, d)

# bring in stage data
stages = read.csv("data/ages2stages.csv")
demog = dplyr::left_join(demog, stages)

demog[which(demog$animal %in% c("L095","L098","L112") & demog$alive==0),] = NA

# write out table of 2021 animals
indx = which(demog$year==2021 & demog$alive==1)

current = demog[indx,]
current = dplyr::select(current, 
                        "animal", "pod", "mom", 
                        "matriline", "sexF1M2", "age","Froh_min10Mb")

current = dplyr::rename(current, froh = Froh_min10Mb)
saveRDS(current, "data/current_pop_2021.rds")

# fit survival models 
demog = dplyr::filter(demog, includeSurv==1,
                    !is.na(alive), birth>1960)
demog$sexF1M2[which(demog$sexF1M2==0)] = sample(1:2,length(which(demog$sexF1M2==0)),replace=T)

# bring in deletrious alleles
del = read.table("data/deleteriousHomozygousCount_srkw", header = TRUE)
names(del) = c("id","homs")
demog = dplyr::left_join(demog, del)