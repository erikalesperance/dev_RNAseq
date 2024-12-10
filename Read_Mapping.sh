#!/bin/bash

# read in the star module
ml star

set -o pipefail

# arrayID (minimum for slurm is zero)
if [[ "${QUEUE}" == "Slurm" ]]; then
	PBS_ARRAYID=$((SLURM_ARRAY_TASK_ID+1))

        # write Job IDs to a text file that will be used to keep track of exit codes
        echo "${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" >> ${RM_JOB_LOG}
	echo "Processing array ${SLURM_ARRAY_TASK_ID} corresponding to line/file # ${PBS_ARRAYID}"
else
	echo "Processing array ${PBS_ARRAYID} through PBS queuing system"
fi


# Get input file (from directory or list)
if [[ -d "$RM_INPUT" ]]; then #if input is a DIRECTORY
	f1=$(find $RM_INPUT -maxdepth 1 -name "*$FORWARD" | sed -n ${PBS_ARRAYID}p)
elif [[ -f "$RM_INPUT" ]]; then #if input is a FILE
	f1=$(sed -n ${PBS_ARRAYID}p $RM_INPUT)
else
	echo "Please specify a valid directory or list of files in the config"
fi

name=$(basename ${f1%%$FORWARD}"") #Name of forward sample

### Get reverse read if paired-end, make sure files are valid
if [[ "$PE" == "True" ]]; then
	f2=${f1%%$FORWARD}"$REVERSE"
	if [[ -f $f1 && -f $f2 ]]; then
		echo "Mapping PE reads for sample $name"
	else
		echo "$f1 and $f2 are not both valid files"
		exit 1
	fi
else
	f2=""
	if [[ -f $f1 ]]; then
		echo "Mapping SE reads for sample $name"
	else
		echo "$f1 is not a valid file"
		exit 1
	fi
fi

### Obtain Sample Name + Lane # to add read group information to mapped files
LANE_NUM=$(grep -o "L00[1-4]" <<< $name | cut -c 4) #Obtain Lane #
SAMPLE_NAME=${name%%[!0-9]*}
ID="${SAMPLE_NAME}:${FLOWCELL_NAME}.${LANE_NUM}"
echo "File name indicates the sample name is ${SAMPLE_NAME} and the lane number is ${LANE_NUM}"
echo "The read group ID field will be ${ID}"

###Genomic Coordinate Output
if [[ "$GENOMIC_COORDINATE_BAMSORTED" == "yes" ]]; then #If genomic alignments should be output as sorted BAM files
	FORMAT="BAM SortedByCoordinate"
	echo "Output Genomic Alignments will be sorted BAM files"
else
	FORMAT="SAM"
	echo "Output Genomic Alignments will be unsorted SAM files"
fi

if [[ "${QUEUE}" == "Slurm" ]]; then
    STAR="STAR"
elif [[ "${QUEUE}" == "PBS" ]]; then
    STAR="${STAR_FILE}"
fi

###star mapping
if [[ "$RM_PASS" == "first" ]]; then ###first pass mode
	echo "In first pass Mode"
	STAR\
	--runThreadN $RM_NTHREAD \
	--genomeDir $GEN_DIR \
	--readFilesIn $f1 $f2 \
	--readFilesCommand gunzip -c \
	--seedSearchStartLmax $SEEDSEARCH \
	--outFileNamePrefix $CJ_OUTPUTDIR/"$name" \
	--outFilterMismatchNmax $MAX_MIS \
	--outFilterMultimapNmax $MAX_N \
	--outFilterScoreMinOverLread $MINSCORE_READL \
	--outFilterMatchNminOverLread $MINMATCH_READL \
	--outReadsUnmapped $UNMAP_F \
	--outSAMtype SAM \
	--quantMode - \
	--outSAMattrRGline ID:${ID} LB:${SAMPLE_NAME} PL:${PLATFORM} SM:${SAMPLE_NAME} PU:${ID} \
	--outFilterType BySJout \
	--outSJfilterReads Unique ## could change later to be a filtering step
elif [[ "$RM_PASS" == "second" ]]; then ###second pass mode
	if [[ ! -z "$JUNCTIONS" ]]; then ### if using junctions
		echo "In second pass mode using $NUM_JUNCTIONS junction files"
		echo "Junctions are as follows: $JUNCTIONS"
		STAR\
		--runThreadN $RM_NTHREAD \
		--genomeDir $GEN_DIR \
		--readFilesIn $f1 $f2 \
		--readFilesCommand gunzip -c \
		--seedSearchStartLmax $SEEDSEARCH \
		--outFileNamePrefix $RM_OUTPUTDIR/"$name" \
		--outFilterMismatchNmax $MAX_MIS \
		--outFilterMultimapNmax $MAX_N \
		--outFilterScoreMinOverLread $MINSCORE_READL \
		--outFilterMatchNminOverLread $MINMATCH_READL \
		--outReadsUnmapped $UNMAP_F \
		--outSAMtype $FORMAT \
		--outSAMattributes NH HI AS nM NM MD \
		--quantMode $QUANT \
		--outSAMattrRGline ID:${ID} LB:${SAMPLE_NAME} PL:${PLATFORM} SM:${SAMPLE_NAME} PU:${ID} \
		--outFilterType BySJout \
		--sjdbFileChrStartEnd $JUNCTIONS
	else
		echo "Mapping without incorporating un-annotated junctions"
		STAR\
		--runThreadN $RM_NTHREAD \
		--genomeDir $GEN_DIR \
		--readFilesIn $f1 $f2 \
		--readFilesCommand gunzip -c \
		--seedSearchStartLmax $SEEDSEARCH \
		--outFileNamePrefix $RM_OUTPUTDIR/"$name" \
		--outFilterMismatchNmax $MAX_MIS \
		--outFilterMultimapNmax $MAX_N \
		--outFilterScoreMinOverLread $MINSCORE_READL \
		--outFilterMatchNminOverLread $MINMATCH_READL \
		--outReadsUnmapped $UNMAP_F \
		--outSAMtype $FORMAT \
		--outSAMattributes NH HI AS nM NM MD \
		--quantMode $QUANT \
		--outSAMattrRGline ID:${ID} LB:${SAMPLE_NAME} PL:${PLATFORM} SM:${SAMPLE_NAME} PU:${ID} \
		--outFilterType BySJout
	fi
else
	echo "Error: Unsure of whether first or second pass mode, exiting..."
	exit 1
fi
