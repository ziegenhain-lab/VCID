#!/bin/bash
# Usage:

# ./run_chunks.sh input.bam input.vcf /path/to/output_dir

# # Requirements:

# 1. samtools, bcftools, GNU parallel, and tabix must be installed (original VCF requires indexing).

# Script Description:

# 1. Use `samtools idxstats` to obtain all chromosome names and lengths, then save them to `chrom_lengths.txt` 
    # in the output path.

# 2. Divide each chromosome into segments of 5,000,000 bases. If the last segment is less than 5,000,000 bases, 
    # use the actual length.

# 3. When extracting from the VCF file, extend the right-hand region by 2,000 bp, but retain the original range 
    # identifiers (e.g., "1-500000") for file and folder names. (Since read pairs are not needed at this stage, 
    # 2,000 bp is sufficient.)

# 4. Only process records with chromosome names of 1-23, X, Y (supporting formats such as "chr1" and "chrX"); 
    # skip other chromosomes.

# 5. Create a separate folder for each segmented region in the specified output path and output BAM and VCF files.

# 6. Use GNU parallel. Execute tasks corresponding to each region in parallel.

# # If you don't need to call Rscript for the time being, you can add a comment before the Rscript command line 
    # (example provided).

# -m means monitor, enabling job monitoring mode, which is Bash's job control mechanism.

# When enabled, Bash can better manage background tasks and process groups.
set -m

BAM_FILE="$1"
VCF_FILE="$2"
OUTPUT_DIR="$3"
READ_SEPARATION_SCRIPT="$4"

if [ -z "$BAM_FILE" ] || [ -z "$VCF_FILE" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$READ_SEPARATION_SCRIPT" ]; then
    echo "Usage: $0 input.bam input.vcf /path/to/output_dir READ_SEPARATION_SCRIPT.R"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Length of each chunk (unit: bp)
CHUNK_SIZE=15000000
# VCF extract region right-side extension bp quantity
VCF_EXTEND=2000
# Please adjust the number of tasks running in parallel according to the actual situation.
NPROC=25

start_time=$(date +%s)

# Use samtools idxstats to obtain chromosomes and their lengths (filtering out lines that start with '*').
samtools idxstats "$BAM_FILE" | awk '$1!="*" {print $1"\t"$2}' > "${OUTPUT_DIR}/chrom_lengths.txt"

# A temporary command file is used to store the commands for all subtasks.
CMD_FILE_MKDIR="${OUTPUT_DIR}/chunk_commands_mkdir.sh"
> "$CMD_FILE_MKDIR"

CMD_FILE_SAM="${OUTPUT_DIR}/chunk_commands_samtools.sh"
> "$CMD_FILE_SAM"

CMD_FILE_VCF="${OUTPUT_DIR}/chunk_commands_bcftools.sh"
> "$CMD_FILE_VCF"

CMD_FILE_READSEPARATION="${OUTPUT_DIR}/chunk_commands_readseparation.sh"
> "$CMD_FILE_READSEPARATION"

# Define a function to determine if a chromosome is valid 
    # (i.e., it conforms to the range 1-23, X, Y, regardless of whether it has the chr prefix).
is_valid_chrom() {
    # Remove any possible "chr" prefixes.
    CHROM_NAME=${1#chr}
    if [[ "$CHROM_NAME" =~ ^(1?[0-9]|2[0-3]|[1-9]|X|Y)$ ]]; then
        return 0
    else
        return 1
    fi
}

cleanup() {
    echo "A termination signal has been detected, and all subprocesses are being killed...."
    kill 0
    exit 1
}

trap cleanup SIGTERM SIGINT


# For each chromosome, perform region segmentation and generate corresponding commands 
    # (only process valid chromosomes, and the output directory is also under the specified output path).
while read -r CHROM LENGTH; do
    if ! is_valid_chrom "$CHROM"; then
        echo "Skip non-specified chromosomes：$CHROM" 
        continue
    fi

    START=1
    while [ "$START" -le "$LENGTH" ]; do
        END=$((START + CHUNK_SIZE - 1))
        # If the current chunk exceeds the total length of the chromosome, then the actual length is used.
        if [ "$END" -gt "$LENGTH" ]; then
            END="$LENGTH"
        fi

        # Calculate the bidirectional extended region
        VCF_START=$((START - VCF_EXTEND))
        # Ensure the starting position is not less than 1.
        if [ "$VCF_START" -lt 1 ]; then
            VCF_START=1
        fi

        # For the VCF extraction region, the right side extends VCF_EXTEND bp.
        VCF_END=$((END + VCF_EXTEND))

        # Ensure the termination position does not exceed the total length of the chromosome.
        if [ "$VCF_END" -gt "$LENGTH" ]; then
            VCF_END="$LENGTH"
        fi

        # Define the directory name of the current chunk, in the format {chromosome}-{start}-{end}.
        CHUNK_DIR="${OUTPUT_DIR}/${CHROM}-${START}-${END}"
        # Define the extracted BAM and VCF file names.
        CHUNK_BAM="${CHROM}_${START}-${END}.bam"
        CHUNK_VCF="${CHROM}_${START}-${END}.vcf"

        # Write the commands for this task into the command file.
        echo "mkdir -p ${CHUNK_DIR}" >> "$CMD_FILE_MKDIR"
        echo "samtools view -b ${BAM_FILE} ${CHROM}:${START}-${END} > ${CHUNK_DIR}/${CHUNK_BAM}" >> "$CMD_FILE_SAM"
        echo "bcftools view -Ov -r ${CHROM}:${VCF_START}-${VCF_END} ${VCF_FILE} > ${CHUNK_DIR}/${CHUNK_VCF}" >> "$CMD_FILE_VCF"
        echo "if [ -f "${CHUNK_DIR}/${CHUNK_BAM}" ] && [ -f "${CHUNK_DIR}/${CHUNK_VCF}" ]; then \
            Rscript ${READ_SEPARATION_SCRIPT} \
                -b ${CHUNK_DIR}/${CHUNK_BAM} \
                -s ${CHUNK_DIR}/${CHUNK_VCF}; \
            fi" >> "$CMD_FILE_READSEPARATION"

        # Update the starting position of the next interval
        START=$((END + 1))
    done
done < "${OUTPUT_DIR}/chrom_lengths.txt"

# Use GNU parallel to execute all tasks in parallel.
# Adding `halt` kills all child processes — when 10 child process errors occur, all child processes are immediately stopped.
parallel --halt now,fail=50 --jobs ${NPROC} < "$CMD_FILE_MKDIR"
parallel --halt now,fail=50 --jobs ${NPROC} < "$CMD_FILE_SAM"
parallel --halt now,fail=50 --jobs ${NPROC} < "$CMD_FILE_VCF"
parallel --halt now,fail=50 --jobs ${NPROC} < "$CMD_FILE_READSEPARATION"

end_time=$(date +%s)
elapsed_seconds=$(( end_time - start_time ))
elapsed_hours=$(echo "scale=2; $elapsed_seconds/3600" | bc)
echo "Total runtime: ${elapsed_hours} hours"
