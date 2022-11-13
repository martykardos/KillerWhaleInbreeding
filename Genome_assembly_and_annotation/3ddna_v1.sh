#!/bin/bash

##########
#The MIT License (MIT)
#
# Copyright (c) 2017 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##########
#
# 3D DNA de novo genome assembly pipeline: version 170123.
#


USAGE_short="
*****************************************************
3D de novo assembly: version 170123

USAGE: ./run-asm-pipeline.sh [options] <path_to_input_fasta> <path_to_input_mnd> <memory (ex:50g)> 

DESCRIPTION:
This is a script to assemble draft assemblies (represented in input by draft fasta and deduplicated list of alignments of Hi-C reads to this fasta as produced by the Juicer pipeline) into chromosome-length scaffolds. The script will produce an output fasta file, Juicebox-compabile contact map of various assembly stages, and a number of supplementary annotation files to help review the result in Juicebox.

ARGUMENTS:
path_to_input_fasta			Specify file path to draft assembly fasta file.
path_to_input_mnd			Specify path to deduplicated list of alignments of Hi-C reads to the draft assembly fasta as produced by the Juicer pipeline: the merged_nodups file (mnd).

OPTIONS:
-m haploid/diploid			Runs in specific mode, either haploid or diploid (default is haploid for Hs2-HiC assembly).
-t tiny_size				Specifies threshold input contig/scaffold size (default is 15000 for Hs2-HiC assembly).
-s number_of_steps			Specifies number of iterative steps for misjoin correction. 
-c chrom_number				Specifies the number of chromosomes expected (default is 23 for for Hs2-HiC assembly).
-j THREAD				Thread for parallel program
-a start_step				The start step for iteration 
-r assembly_mode			The mode for assembly,[scaffold/contig/boundary]
-S stage                                Fast forward to later assembly steps, can be init,iter,split,seal,final,agp,tohic,heatmap,sca_len,n50,report 
-h					Shows this help. Type -e for a longer set of options.

I'm so sorry that some parameters have no relative meaning with their role,because there are so many parameters.


*****************************************************
"

USAGE_long="
*****************************************************
3D de novo assembly: version 170123

USAGE: ./run-asm-pipeline.sh [options] <path_to_input_fasta> <path_to_input_mnd> 

DESCRIPTION:
This is a script to assemble draft assemblies (represented in input by draft fasta and deduplicated list of alignments of Hi-C reads to this fasta as produced by the Juicer pipeline) into chromosome-length scaffolds. The script will produce an output fasta file, Juicebox-compabile contact map of various assembly stages, and a number of supplementary annotation files to help review the result in Juicebox.

ARGUMENTS:
path_to_input_fasta			Specify file path to draft assembly fasta file.
path_to_input_mnd			Specify path to deduplicated list of alignments of Hi-C reads to the draft assembly fasta as produced by the Juicer pipeline: the merged_nodups file (mnd).

OPTIONS:
-e					Shows this help.
-m haploid/diploid			Runs in specific mode, either haploid or diploid (default is haploid). In this version of the pipeline diploid assumes Rabl.
-s number_of_steps			Specifies number of iterative steps for misassembly detector (default is 2).
-t tiny_size				Specifies threshold input contig/scaffold size for initial filtration (default is 15000 or 15kb). We do not attempt to further assembly contigs/scaffolds smaller than this size.
-c chrom_number				Number of chromosomes expected (default is 23 for human assembly).

OPTIONS (hardcoded at the moment):
**scaffolder**
-q mapq					Mapq threshold for scaffolding (default is 1).

**misjoin detector**
-w r1_bin_size				Misjoin detection coarse matrix resolution, should be one of the following: 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000 (default is 25000 or 25kb).
-n r2_bin_size				Misjoin detection fine matrix resolution, should be one of the following: 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000 (default is 1000 or 1kb).
-d d_region_size			Misjoin detection triangular motif region size (default is 125000 or 125kb).
-k stringency				Misjoin detection stringency parameter (default is 55).
-p percentile				Misjoin detection saturation parameter (default is 5).

**polisher**
-x r1_bin_size				Large-scale misjoin detection coarse matrix resolution for Rabl polishing (default is 100000 or 100kb).
-y d_region_size			Large-scale misjoin detection coarse matrix resolution for Rabl polishing (default is 3000000 or 3Mb).
-z stringency				Large-scale misjoin detection stringency parameber for Rabl polishing (default is 55).

**merger**
-b band_size				Band size for alternative haplotype detection (default is 1000000 or 1Mb).
-l lastz_opt_string			LASTZ aligner option string.
-i identity					Sequence identity threshold.

**workflow**
-S stage					Fast forward to later assembly steps, can be polish, split, seal, merge and finalize
*****************************************************
"

pipeline=$(dirname $(readlink -f $0))
#source $pipeline/config_path.txt
##++++++++++ default parameter setup+++++++++++++++
diploid="false"	# by default run haploid pipeline
input_size=15000	# tiny contig size threshold; used 15kb for Hs2 and 10kb for AaegL4
mapq=1	# read mapping quality threshold for Hi-C scaffolder
THREAD=1
# misassembly detection default params
wide_bin=25000	
narrow_bin=1000
wide_depletion=125000
sensitivity=55
pct=5
thr_cov=2

############### HANDLE OPTIONS ###############

while getopts "hem:t:s:c:a:q:w:n:d:k:p:S:j:r:" opt; do
case $opt in
	h) echo "$USAGE_short" >&1
		exit 0
	;;
	e) echo "$USAGE_long" >&1
		exit 0
	;;
## short menu
	m)	if [ "$OPTARG" == "haploid" ] || [ "$OPTARG" == "diploid" ]; then
    	    echo >&1 "-m flag was triggered. Running in $OPTARG mode."
			[ "$OPTARG" == "diploid" ] && diploid="true" && thr_cov=4	# use higher repeat coverage for tentatively heterozygous sequence
    	else
    		echo ":( Unrecognized value for -m flag. Running with default parameters (-m haploid)." >&2
    	fi
    ;;
	t)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -t flag was triggered, filtering draft contigs/scaffolds smaller than $OPTARG" >&1
            input_size=$OPTARG
        else
            echo ":( Wrong syntax for tiny size threshold. Using the default value ${input_size}" >&2
        fi
    ;;
    s)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -ge 0 ]]; then
            echo ":) -s flag was triggered, will run $OPTARG rounds of misjoin correction" >&1
            MAX_STEP=$OPTARG
        else
            echo ":( Wrong syntax for number of iterative steps for misjoin correction. Using the default value ${MAX_STEP}" >&2
        fi
	;;
	c)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]]; then
            echo ":) -c flag was triggered, assuming $OPTARG chromosomal-scaffolds in output" >&1
            chrom_num=$OPTARG
        else
            echo ":( Wrong syntax for chromosome number. Using the default value ${chrom_num}" >&2
        fi
	;;
	a) re='^[0-9]+$'
	if [[ $OPTARG =~ $re ]]; then
	    start_step=$OPTARG
    	fi
	;;	
## extended menu: scaffolder
	q)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -q flag was triggered, scaffolding using reads with at least $OPTARG mapping quality" >&1
            mapq=$OPTARG
        else
            echo ":( Wrong syntax for mapping quality. Using the default value ${mapq}" >&2
        fi
    ;;
## extended menu: misjoin detector
    w)  re='^[0-9]+$'	## TODO: specify/generalize re matrix resolutions size
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -w flag was triggered, misjoin detection coarse matrix resolution set to $OPTARG" >&1
            wide_bin=$OPTARG
        else
            echo ":( Wrong syntax for misjoin detection coarse matrix resolution. Using the default value ${wide_bin}" >&2
        fi
    ;;
    n)  re='^[0-9]+$'	## TODO: specify/generalize re matrix resolutions size
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -n flag was triggered, misjoin detection fine matrix resolution set to $OPTARG" >&1
            narrow_bin=$OPTARG
        else
            echo ":( Wrong syntax for misjoin detection fine matrix resolution. Using the default value ${narrow_bin}" >&2
        fi
    ;;
    d)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -d flag was triggered, misjoin detection depletion region size set ot $OPTARG" >&1
            wide_depletion=$OPTARG
        else
            echo ":( Wrong syntax for misjoin detection depletion region size. Using the default value ${wide_depletion}" >&2
        fi
    ;;
    k)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]]; then
            echo ":) -k flag was triggered, misjoin detection stringency parameter set ot $OPTARG" >&1
            sensitivity=$OPTARG
        else
            echo ":( Wrong syntax for misjoin detection stringency parameter. Using the default value ${sensitivity}" >&2
        fi
    ;;
    p)  re='^[0-9]+$'
        if [[ $OPTARG =~ $re ]] && [[ $OPTARG -gt 0 ]] && [[ $OPTARG -lt 100 ]]; then
            echo ":) -p flag was triggered, misjoin detection saturation parameter set to ${OPTARG}%" >&1
            pct=$OPTARG
        else
            echo ":( Wrong syntax for misjoin detection saturation threshold. Using the default value ${pct}%" >&2
        fi
    ;;
	S) #if [ "$OPTARG" == "scaffold" ] || [ "$OPTARG" == "polish" ] || [ "$OPTARG" == "split" ] || [ "$OPTARG" == "seal" ] || [ "$OPTARG" == "merge" ] || [ "$OPTARG" == "finalize" ]; then
	#		echo ":) -S flag was triggered, fast-forwarding to \"$OPTARG\" pipeline section" >&1
			stage=$OPTARG
	#	else
#			echo ":( Wrong syntax for pipeline stage. Exiting!" >&2
#		fi		
	;;
        j) re='^[0-9]+$'
            if [[ $OPTARG =~ $re ]]; then
                echo ":) -j flag was triggered" >&1
                THREAD=$OPTARG
            else
                echo ":( Wrong syntax for threads" >&2
            fi
        ;;
        r) echo ":) -r flag was triggered" >&1		
            assembly_mode=$OPTARG
	;;    
	*)  echo ":( Illegal options. Exiting."
		echo "$USAGE_short"
		exit 1
	;;
esac
done

shift $(( OPTIND-1 ))

##++++++++++++++++++++++check parameters for compatibility+++++++++++++++++++
[[ ${wide_depletion} -le ${wide_bin} ]] && echo >&2 ":( Requested depletion region size ${wide_depletion} and bin size ${wide_bin} parameters are incompatible. Run ${pipeline}/edit/run-mismatch-detector.sh -h for instructions. Exiting!" && exit 1

[[ ${wide_bin} -le ${narrow_bin} ]] && echo >&2 ":( Requested mismatch localization resolution ${narrow_bin} and coarse search bin size ${wide_bin} parameters are incompatible. Run ${pipeline}/edit/run-mismatch-detector.sh -h for instructions. Exiting!" && exit 1

##++++++++++++++++++++++check parameters for consistency and set up defaults for diploid mode++++++++++++++++
if [ "$diploid" == "true" ]; then
	[ -z ${polish_wide_bin} ] && polish_wide_bin=100000
	[ -z ${polish_wide_depletion} ] && polish_wide_depletion=3000000
	[ -z ${polish_sensitivity} ] && polish_sensitivity=55
	[ -z ${polish_size} ] && polish_size=1000000
	[ -z ${band_size} ] && band_size=1000000
fi

#++++++++++++++++++++ HANDLE EXTERNAL DEPENDENCIES +++++++++++++++++++++++++ 

##	GNU Parallel Dependency
parallel="false"
if hash parallel 2>/dev/null; then
        ver=`parallel --version | awk 'NR==1{print \$3}'`
        [ $ver -ge 20150322 ] && parallel="true"
fi
[ $parallel == "false" ] && echo ":| WARNING: GNU Parallel version 20150322 or later not installed. We highly recommend to install it to increase performance. Starting pipeline without parallelization!" >&2

## LASTZ dependency
lastz="false"
if hash lastz 2>/dev/null; then
	lastz="true"
fi
[[ $lastz == "false" && $diploid == "true" ]] && echo >&2 ":( LASTZ not installed or not in the path while diploid mode is triggered. The merge section of the pipeline will not run. Exiting!" && exit 1

#++++++++++++++++++ HANDLE ARGUMENTS +++++++++++++++++++++++++

orig_fasta=$1
orig_mnd=$2
memory=$3
genomeid=$(basename "$orig_fasta" .fasta)
genomeid=$(basename "$genomeid" .fna)
genomeid=$(basename "$genomeid" .fa)

#+++++++++++++init stage++++++++++++++++++++
if [[ $diploid == "false" ]];then
    stage_all="init,iter,split,seal,final,agp,tohic,heatmap,sca_len,n50,report"
else
    stage_all="init,iter,polish,split,seal,merge,final,agp,tohic,heatmap,sca_len,n50,report"
fi

#+++++++++++++stage for continue++++++++++++++++++

if [[ $stage ]];then
    if [[ $stage =~ "manu_split" ]];then
        stage_all=$(echo $stage_all|sed 's/split/manu_split/g')
    fi
    if [[ $stage =~ "-" ]];then
        stage=$(echo $stage|sed 's/-//g')
        stage_all=$(echo $stage_all|sed s/.*$stage//g)
        stage_all=$stage""$stage_all
    else
        stage_all=$stage
    fi
fi
echo $stage_all
#+++++++++++++++++++init++++++++++++++++++++++++++
orig_cprops=${genomeid}.cprops
if [[ $stage_all =~ "init" ]];then
    [ -z $1 ] || [ -z $2 ] && echo >&2 "Not sure how to parse your input: files not listed or not found at expected locations. Exiting!" && echo >&2 "$USAGE_short" && exit 1

    [ ! -s $1 ] || [ ! -s $2 ] && echo >&2 "Not sure how to parse your input: files not listed or not found at expected locations. Exiting!" && echo >&2 "$USAGE_short" && exit 1

## TODO: check file format
    if [ "$#" -ne 3 ]; then
        echo >&2 "Illegal number of arguments. Please double check your input. Exiting!" && echo >&2 "$USAGE_short" && exit 1
    fi
    #+++++++++++generate mnd++++++++++++++
    ln -sf ${orig_mnd} ${genomeid}.mnd.txt
    orig_mnd=${genomeid}.mnd.txt
    
    #+++++++++++++generate cprops+++++++++++++++
    awk -f ${pipeline}/utils/generate-cprops-file.awk ${orig_fasta} > ${genomeid}.cprops
#    orig_cprops=${genomeid}.cprops
    
    #+++++++++++++++check cprops++++++++++++++++
    [ ! -f ${orig_cprops} ] && echo >&2 ":( No cprops file found. Please rerun the pipeline from scratch. Exiting!" && exit 1
    #+++++++++++++++calculate scale+++++++++++++++
fi
if [ -f ${orig_cprops} ];then
    totlength=`awk '{total+=$3}END{print total}' ${orig_cprops}`
    scale=$(( 1 + $totlength / 2100000000 ))
    if [[ $scale -ne 1 ]]; then
	    wide_bin=$((wide_bin/scale))
	    wide_depletion=$((wide_depletion/scale))
	    narrow_bin=$((narrow_bin/scale))
    fi
    ln -sf ${orig_cprops} ${genomeid}.0.cprops
    ln -sf ${orig_mnd} ${genomeid}.mnd.0.txt
fi
#ln -sf ${orig_cprops} ${genomeid}.0.cprops
#ln -sf ${orig_mnd} ${genomeid}.mnd.0.txt

#++++++++++++++++++++++ ITERATIVE SCAFFOLDING/MISJOIN CORRECTION +++++++++++++++++
if [[ $stage_all =~ "iter" ]];then
    
#	touch ${genomeid}.edits.txt

    if [[ $assembly_mode && $assembly_mode == "scaffold" ]];then
        start_step=0
        MAX_STEP=0
    elif [[ $assembly_mode && $assembly_mode == "contig" ]];then
        python ${pipeline}/scripts/find_N.py ${orig_fasta} > N_loc.txt
    fi
    
    if [[ $start_step -eq 0 ]];then
        touch "archive."${genomeid}".edits.at.step.0.txt"
    fi
   
    if [[ ! $start_step ]];then
        start_step=0
    fi
    current_cprops=${genomeid}.${start_step}.cprops
    current_mnd=${genomeid}.mnd.${start_step}.txt

#    STEP=$start_step	
    if [[ ! -f "archive."${genomeid}".edits.at.step."$((start_step))".txt" ]];then
        echo "Error:missing archive."${genomeid}".edits.at.step."$((start_step))".txt,maybe the step was setted error!" && exit 1
    else
        cp "archive."${genomeid}".edits.at.step."$((start_step))".txt" ${genomeid}".edits.txt"
    fi
    STEP=$start_step
    while true; do
        echo $STEP 
        if [[ $start_step -eq 0 || $start_step -lt $STEP ]];then

		    bash ${pipeline}/scaffold/run-liger-scaffolder.sh -p ${parallel} -j ${THREAD} -s ${input_size} -q ${mapq} ${current_cprops} ${current_mnd}
        	# build a hic map of the resulting assembly
		    bash ${pipeline}/visualize/run-asm-visualizer.sh -p ${parallel} -j ${THREAD} -q ${mapq} ${current_cprops} ${genomeid}.${STEP}.asm ${current_mnd} ${STEP} $((MAX_STEP+1)) ${memory}
		    rm temp.${genomeid}.${STEP}.asm_mnd.txt

        fi
	# break out of the scaffold-mismatch detection loop if the max number of steps is reached
		[ ${STEP} -eq ${MAX_STEP} ] && break

	# annotate near-diagonal mismatches in the map
		bash ${pipeline}/edit/run-mismatch-detector.sh -p ${parallel} -j ${THREAD} -w ${wide_bin} -d ${wide_depletion} -n ${narrow_bin} -k ${sensitivity} ${genomeid}.${STEP}.hic

	# annotate repeats by coverage analysis
		bash ${pipeline}/edit/run-coverage-analyzer.sh -w ${wide_bin} -t ${thr_cov} ${genomeid}.${STEP}.hic

	# store intermediate mismatch stuff	- not necessary
		mv depletion_score_wide.wig depletion_score_wide.at.step.${STEP}.wig
		mv depletion_score_narrow.wig depletion_score_narrow.at.step.${STEP}.wig
		mv mismatch_wide.bed mismatch_wide.at.step.${STEP}.bed
		mv mismatch_narrow.bed mismatch_narrow.at.step.${STEP}.bed
	# store intermediate repeat stuff - not necessary
		mv coverage_wide.wig coverage_wide.at.step.${STEP}.wig
		mv repeats_wide.bed repeats_wide.at.step.${STEP}.bed

	# consolidate bed annotations
		cat mismatch_narrow.at.step.${STEP}.bed repeats_wide.at.step.${STEP}.bed | sort -k 2,2n | awk 'BEGIN{FS="\t"; OFS="\t"}NR==1{start=$2; end=$3; next}$2<=end{if($3>end){end=$3}; next}{print "assembly", start, end; start=$2; end=$3}END{print "assembly", start, end}' > suspect.at.step.${STEP}.bed

	# convert bed track into 2D annotations
		resolved=$(awk 'NR==2{print $3}' ${genomeid}.${STEP}_asm.superscaf_track.txt)	# scaled coordinates	
		awk -v bin_size=${narrow_bin} -f ${pipeline}/edit/overlay-edits.awk ${genomeid}.${STEP}_asm.scaffold_track.txt suspect.at.step.${STEP}.bed | awk -v r=${resolved} 'NR==1||$3<=r' > suspect_2D.at.step.${STEP}.txt

	# separate intra and inter-input scaffold mismatches
		awk 'NR==1||$8=="debris"' suspect_2D.at.step.${STEP}.txt > edits.for.step.$((STEP+1)).txt
	# optional
		awk 'NR==1||$8=="mismatch"' suspect_2D.at.step.${STEP}.txt > mismatches.at.step.$STEP.txt
		
		test=`wc -l < edits.for.step.$((STEP+1)).txt`
	
		[ $test -eq 1 ] && echo >&1 ":) No more input edits to be done. Moving to polishing!" && break

	# move on to the next step	
		STEP=$((STEP+1))
		cp ${genomeid}".edits.txt" "archive."${genomeid}".edits.at.step."$((STEP-1))".txt"
	
		bash ${pipeline}/lift/lift-edit-asm-annotations-to-original-input-annotations.sh ${orig_cprops} ${current_cprops} ${genomeid}.$((STEP-1)).asm edits.for.step.${STEP}.txt > h.edits.txt
		awk 'NR==1' "h.edits.txt" > temp
		{ awk 'NR>1' "archive."${genomeid}".edits.at.step."$((STEP-1))".txt" ; awk 'NR>1' "h.edits.txt" ; } | sort -k 1,1 -k 2,2n >> temp
		mv temp ${genomeid}".edits.txt"

		[[ $assembly_mode && $assembly_mode == "contig" ]] && mv ${genomeid}".edits.txt" ${genomeid}".edits.txt.bak" && python ${pipeline}/scripts/edit_edits.py N_loc.txt ${genomeid}".edits.txt.bak" > ${genomeid}".edits.txt"	

		bash ${pipeline}/edit/apply-edits-prep-for-next-round.sh -p ${parallel} -j ${THREAD} -r ${STEP} ${genomeid}".edits.txt" ${orig_cprops} ${orig_mnd}
		current_cprops=${genomeid}.${STEP}.cprops
		current_mnd=${genomeid}.mnd.${STEP}.txt
	
	done

	ln -sf ${genomeid}.${STEP}.cprops ${genomeid}.resolved.cprops
	ln -sf ${genomeid}.${STEP}.asm ${genomeid}.resolved.asm
	ln -sf ${genomeid}.${STEP}_asm.scaffold_track.txt ${genomeid}.resolved_asm.scaffold_track.txt
	ln -sf ${genomeid}.${STEP}_asm.superscaf_track.txt ${genomeid}.resolved_asm.superscaf_track.txt
	ln -sf ${genomeid}.${STEP}.hic ${genomeid}.resolved.hic
	ln -sf ${genomeid}.mnd.${STEP}.txt ${genomeid}.mnd.resolved.txt
	
	if [ $diploid == "false" ]; then
		ln -sf ${genomeid}.resolved.cprops ${genomeid}.polished_pro.cprops
		ln -sf ${genomeid}.resolved.asm ${genomeid}.polished_pro.asm
		ln -sf ${genomeid}.mnd.resolved.txt ${genomeid}.mnd.polished.txt 
	fi

fi

############### POLISHING ###############

if [[ $stage_all =~ "polish" ]];then
	[ ! -f ${genomeid}.resolved.cprops ] || [ ! -f ${genomeid}.resolved.asm ] || [ ! -f ${genomeid}.resolved.hic ] || [ ! -f ${genomeid}.mnd.resolved.txt ] && echo >&2 ":( No resolved files are found. Please rerun the pipeline to include the scaffold segment. Exiting!" && exit 1
	if [ $diploid == "true" ]; then
		bash ${pipeline}/polish/run-asm-polisher.sh -p ${parallel} -j ${THREAD} -g ${genomeid}.resolved.hic -a ${genomeid}.resolved_asm.scaffold_track.txt -b ${genomeid}.resolved_asm.superscaf_track.txt -s ${polish_size} -w ${polish_wide_bin} -d ${polish_wide_depletion} -n ${narrow_bin} -k ${polish_sensitivity} ${genomeid}.cprops ${orig_mnd} ${genomeid}.resolved.cprops ${genomeid}.resolved.asm
		
		rm temp_resolved.asm temp_unresolved.asm new_resolved.asm new_unresolved.asm temp.pre_polish_edits.txt temp.post_polish_edits.txt temp.${genomeid}.resolved.polish.asm_mnd.txt
		
		mv ${genomeid}.resolved.polish.cprops ${genomeid}.polished_pro.cprops
		mv ${genomeid}.resolved.polish.asm ${genomeid}.polished_pro.asm

		mv ${genomeid}.resolved.polish.edits_2D.txt ${genomeid}.polished.edits_2D.txt 
		mv ${genomeid}.resolved.polish.mismatches_2D.txt ${genomeid}.polished.mismatches_2D.txt
		mv ${genomeid}.resolved.polish.suspect_2D.txt ${genomeid}.polished.suspect_2D.txt
		mv ${genomeid}.resolved.polish.mismatch_narrow.bed ${genomeid}.polished.mismatch_narrow.bed
		mv ${genomeid}.resolved.polish.depletion_score_narrow.wig ${genomeid}.polished.depletion_score_narrow.wig
		mv ${genomeid}.resolved.polish.mismatch_wide.bed ${genomeid}.polished.mismatch_wide.bed
		mv ${genomeid}.resolved.polish.depletion_score_wide.wig ${genomeid}.polished.depletion_score_wide.wig
		
		mv ${genomeid}.resolved.polish.hic ${genomeid}.polished.hic
		mv ${genomeid}.resolved.polish_asm.superscaf_track.txt ${genomeid}.polished_asm.superscaf_track.txt
		mv ${genomeid}.resolved.polish_asm.scaffold_track.txt ${genomeid}.polished_asm.scaffold_track.txt
		mv ${genomeid}.resolved.polish.mnd.txt ${genomeid}.mnd.polished.txt

	fi
fi
############### SPLITTING ###############

if [[ $stage_all =~ "split" ]];then
    STEP=${MAX_STEP}
    if [ ! -f "juicebox_assembly.txt" ];then
        [ $chrom_num -ne 1 ] && bash ${pipeline}/split/run-asm-splitter.sh -c ${chrom_num} -r ${diploid} ${genomeid}.polished_pro.cprops ${genomeid}.polished_pro.asm ${genomeid}.mnd.polished.txt || cp ${genomeid}.polished.asm ${genomeid}.polished.split.asm
    fi
    split_pro_num=$(sed 's/ /\n/g' ${genomeid}.polished_pro.split.asm|wc -l) 
    split_num=$(( $split_pro_num ))
    orgin_num=$(sed 's/ /\n/g' ${genomeid}.${STEP}.asm|wc -l)
    if [ ${split_num} -eq ${orgin_num} ];then
        cat <(awk '{$0=">"$0}1' ${genomeid}.polished_pro.cprops) ${genomeid}.polished_pro.split.asm > ${genomeid}.${STEP}.assembly 
    else
        cat <(awk '{$0=">"$0}1' ${genomeid}.${STEP}.cprops) ${genomeid}.${STEP}.asm > ${genomeid}.${STEP}.assembly
    fi

    while :
    do
        if [ -f "juicebox_assembly.txt" ];then
            dos2unix juicebox_assembly.txt
            break
        fi
        sleep 10m
    done

    awk -v cprops=${genomeid}.polished.cprops -v asm=${genomeid}.polished.split_pro.asm '$1~/^>/{$1=substr($1,2); print > cprops;next}{print > asm}' juicebox_assembly.txt
    awk -v chr_num=${chrom_num} 'NR>chr_num{gsub(" ","\n",$0);print $0;next}{print $0}'  ${genomeid}.polished.split_pro.asm >${genomeid}.polished.split.asm
	
fi

#+++++++++++++ manu_split ++++++++++++++++

if [[ $stage_all =~ "manu_split" ]];then
    python ${pipeline}/utils/manu_split.py ${genomeid}.polished.asm ${genomeid}.polished.cprops "boundaries.list" && cp split_asm ${genomeid}.polished.split.asm 
fi

############### SEALING ###############
if [[ $stage_all =~ "seal" ]];then [ ! -s ${genomeid}.polished.cprops ] || [ ! -s ${genomeid}.polished.split.asm ] && echo >&2 ":( No split files are found. Please rerun the pipeline to include the split segment. Exiting!" && exit 1
	
	bash ${pipeline}/seal/seal-asm_v1.sh -s ${input_size} ${orig_cprops} ${genomeid}.polished.cprops ${genomeid}.polished.split.asm ${chrom_num} 
        # change seal-asm.sh to seal-asm_v1.sh and add chrom_num  by qun 20201027

	mv ${genomeid}.polished.sealed.cprops ${genomeid}.rawchrom.cprops
	mv ${genomeid}.polished.sealed.asm ${genomeid}.rawchrom.asm

	## TODO: build a hic map here: downside takes some time.

	# prep for merging and finalizing
	python ${pipeline}/utils/Toraw_fasta.py ${genomeid}.rawchrom.cprops ${orig_fasta} ${genomeid}.rawchrom.fasta
	

	
	if [ $diploid == "false" ]; then
		ln -sf ${genomeid}.rawchrom.cprops ${genomeid}.final.cprops
		ln -sf ${genomeid}.rawchrom.asm ${genomeid}.final.asm
		ln -sf ${genomeid}.rawchrom.fasta ${genomeid}.final.fasta
	fi
	
fi

############### MERGING ###############
if [[ $stage_all =~ "merge" ]];then
	
	[ ! -s ${genomeid}.rawchrom.cprops ] || [ ! -s ${genomeid}.rawchrom.asm ] || [ ! -s ${genomeid}.rawchrom.fasta ] && echo >&2 ":( No raw chromosomal files were found. Please rerun he pipeline to include the seal segment" && exit 1
	
	echo ":) Starting merge block" >&1
	
## TODO: split unsafe, redo via indexing as in haploid case
	[ -d faSplit ] && echo >&2 "faSplit directory exists. Remove faSplit directory and restart with -S merge flag. Exiting!" && exit 1
	mkdir faSplit && cd faSplit && awk -f ${pipeline}/merge/split-fasta-by-cname.awk ../${genomeid}.rawchrom.cprops ../${genomeid}.rawchrom.fasta && cd ..
		
	bash ${pipeline}/merge/run-asm-merger.sh -b ${band_size} ${genomeid}.rawchrom.cprops ${genomeid}.rawchrom.asm faSplit
	[ -f joblist.txt ] && rm joblist.txt

	cp ${genomeid}.rawchrom/${genomeid}.rawchrom_merged.asm ${genomeid}.final.asm
	mv ${genomeid}.rawchrom/merged_${genomeid}.rawchrom.fa ${genomeid}.final.fasta
	awk -f ${pipeline}/utils/generate-cprops-file.awk ${genomeid}.final.fasta > ${genomeid}.final.cprops
	
fi

#+++++++++++++++ final ++++++++++++++
if [[ $stage_all =~ "final" ]];then

    #+++++++++++++++++  clean.sh  +++++++++++++++
    /bin/ls |grep -Ev "*.FINAL.from_draft.fasta|merged_nodups.txt|*\.sh|*\.sh\.e*|*\.sh\.o*|boundary.list|${genomeid}.${STEP}.hic|*.asm|*.cprops"|awk '{print "/bin/rm "$1}' > clean.sh 
    #+++++++++++++++++ resort asm +++++++++++++++
    if [[ ! $stage_all =~ "manu_split" ]];then
        ${pipeline}/utils/resortAsm.pl -l ${genomeid}.final.cprops -i ${genomeid}.final.asm -n ${chrom_num} 
    fi
    #+++++++++++++++++ FINALIZING +++++++++++++++
    bash ${pipeline}/finalize/finalize-output_v3.sh -c ${chrom_num} -s ${input_size} -l ${genomeid} ${genomeid}.final.cprops ${genomeid}.final.asm ${genomeid}.final.fasta final

fi
############### CLEANUP ###############

# only big files, will ultimately remove more intermediate stuff TODO: remove in respective blocks
#rm ${genomeid}.mnd.*.txt
[ -d faSplit ] && rm faSplit
[ -d ${genomeid}.rawchrom ] && rm -r ${genomeid}.rawchrom && rm -r ${genomeid}.final.fasta

#++++++++++++++ agp +++++++++++++++++++++
if [[ $stage_all == "agp" ]];then
    if [[ $# -lt 5 ]];then
        echo "bash /script_path/3d_dna.sh -S agp orig_fasta cprops asm chr_number out\n" && exit 1
    fi
    orig_fasta=$1
    cprops=$2
    asm=$3
    chr_num=$4
    out=$5
elif [[ $stage_all =~ "agp" ]];then
    orig_fasta=${orig_fasta}
    cprops=${genomeid}.final.cprops
    asm=${genomeid}.final.asm
    chr_num=${chrom_num}
    out=${genomeid}_agp.out
fi
if [[ $stage_all =~ "agp" ]];then
    ${pipeline}/utils/3d-dna2agp.v2.5.pl -draft ${orig_fasta} -cprop ${cprops} -asm ${asm} -chr ${chr_num} -out ${out} -tiny_size ${input_size} 
fi
#+++++++++++++ tohic +++++++++++++++++++++
if [[ $stage_all == "tohic" ]];then
    if [[ $# -ne 5 ]];then
        echo "bash /script_path/3d_dna.sh -S tohic cprops origin_mnd(ex:merged_nodup.txt) asm cpu memory(100g)\n" && exit 1
    fi
    cprops=$1
    mnd=$2
    asm=$3
    cpu=$4
    mem=$5
elif [[ $stage_all =~ "tohic" ]];then
    cprops=${genomeid}.final.cprops
    mnd=${orig_mnd}
    asm=${genomeid}.final.asm
    cpu=${THREAD}
    mem=${memory}
fi
if [[ $stage_all =~ "tohic" ]];then
    bash ${pipeline}/utils/final2hic.sh $cprops $mnd $asm $cpu 1 2 $mem
    cat <(awk '{$0=">"$0}1' $cprops) $asm > final.assembly
fi
#+++++++++++++ heatmap ++++++++++++++++++++++
if [[ $stage_all == "heatmap" ]];then
    if [[ $# -ne 4 ]];then
        echo "bash /script_path/3d_dna.sh -S heatmap hic cprops asm chr_num\n" && exit 1
    fi
    hic=$1
    cprops=$2
    asm=$3
    chr_num=$4
elif [[ $stage_all =~ "heatmap" ]];then
    hic=${genomeid}.final.hic
    cprops=${genomeid}.final.cprops
    asm=${genomeid}.final.asm
    chr_num=${chrom_num}
fi
if [[ $stage_all =~ "heatmap" ]];then
    bash ${pipeline}/utils/final2draw.sh $hic $cprops $asm $chr_num 1 1
fi
#+++++++++++ length of scaffold ++++++++++
if [[ $stage_all == "sca_len" ]];then
    if [[ $# -ne 2 ]];then
        echo "bash /script_path/3d_dna.sh -S sca_len genome out"
    fi
    genome=$1
    out=$2
elif [[ $stage_all =~ 'sca_len' ]];then
    genome=${genomeid}.FINAL.from_draft.fasta
    out=scaffold_length.txt
fi
if [[ $stage_all =~ 'sca_len' ]];then
    awk '$1~/^>/{if(seq_len != 0){print name"\t"seq_len};name=$1;seq_len=0;next}{seq_len+=length($1)}END{print name"\t"seq_len}' $genome > $out
fi
#++++++++++++++ n50 ++++++++++++++++++++
if [[ $stage_all == "n50" ]];then
    if [[ $# -ne 2 ]];then
        echo "bash /script_path/3d_dna.sh -S n50 genome out"
    fi
    genome=$1
    out=$2
elif [[ $stage_all =~ 'n50' ]];then
    genome=${genomeid}.FINAL.from_draft.fasta
    out=N50.txt
fi
if [[ $stage_all =~ 'n50' ]];then
    ${pipeline}/utils/cs  ${genome} > ${out}
fi
#+++++++++++++++ report ++++++++++++++++
if [[ $stage_all =~ 'report' ]];then
    mkdir 3ddna_report
    mv N50.txt 3ddna_report/
    mv *final.hic 3ddna_report/
    mv scaffold_length.txt 3ddna_report/
    cd 3ddna_report
    ln -s ../*FINAL.from_draft.fasta .
cd ..
fi
echo "Mission Complete" > job_status.log
