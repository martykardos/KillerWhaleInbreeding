README for folder ROH

Filtering information for the data used in the ROH analysis, and the sequence of the analyses are described in orca_filterSNPs_estimateNe


popKey                                                            a key showing which population different killer whale IDs belong to.
rCode_prepareGenotypeProbabilitiesForROH_1June2022.R              prepares genotype likelihoods for use in identifying ROH. Input is 
rCode_filterLociHardyDepth_22July2022.R                           filters loci for ROH analysis based on read depth and deviation from Hardy Weinberg proportions
rCode_lodRohDetection_genoLikelihoods_2June2022.R                 identifies runs of homozygosity in the genome sequence data
rCode_calculateFroh_5June2022.R                                   calculates F_ROH using ROH identified in 


The sequence of the analysis is outlined in the file "orca_filterSNPs_estimateNe"