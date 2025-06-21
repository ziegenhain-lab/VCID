### read_mapping_and_bam_sorting_sub_pipeline.smk 子工作流 ----------------------------------------
import yaml
import os

# 1.读取传递的 config_file -----------------------------------------------------------
config_file = config["config_file"]  
with open(config_file, 'r') as f:
    batch_config = yaml.safe_load(f)


# 2.Input files from batch config 基础信息 ----------------------------------------------
PROJECT_NAME = batch_config['project']
Program_folder = batch_config['pipeline_path']
algorithem_path = os.path.join(Program_folder, "algorithem", "snakemake_pipeline")
CELLTYPE1 = batch_config["celltype1"]
CELLTYPE2 = batch_config["celltype2"]
CELLTYPE3 = batch_config["celltype3"]
ANNOTATION = batch_config['annotation']

Data_results_path = batch_config['out_dir']
OUT_DIR = Data_results_path
os.makedirs(OUT_DIR, exist_ok=True)

GTF_FILE = batch_config['reference']['GTF_file']
REFERENCE_GENOME = batch_config['reference']['reference_genome']
REFERENCE_GENOME_INDEX = f"{REFERENCE_GENOME}.fai"
THREADS = batch_config['num_threads']
ZUMIS_PATH = batch_config['zUMIs_directory']

# 3. Output files  ----------------------------------------------------------------
ZUMIS_OUTPUT_DIR = os.path.join(OUT_DIR, "zUMIs_output")
os.makedirs(ZUMIS_OUTPUT_DIR, exist_ok=True)

SPLIT_BAM_OUTPUT_DIR = os.path.join(OUT_DIR, "split_bam_output")
os.makedirs(SPLIT_BAM_OUTPUT_DIR, exist_ok=True)

READ_SEPARATION_RESULTS_DIR = os.path.join(OUT_DIR, "read_separation_results")
os.makedirs(READ_SEPARATION_RESULTS_DIR, exist_ok=True)


# 4.Rule生成的文件的信息标注  -----------------------------------------------------
# 1）Extract celltype specific barcode 产生的文件
CELLBC_CELLTYPE1 = os.path.join(SPLIT_BAM_OUTPUT_DIR, f"{CELLTYPE1}_subset_cellBClist.txt")
CELLBC_CELLTYPE2 = os.path.join(SPLIT_BAM_OUTPUT_DIR, f"{CELLTYPE2}_subset_cellBClist.txt")
CELLBC_CELLTYPE3 = os.path.join(SPLIT_BAM_OUTPUT_DIR, f"{CELLTYPE3}_subset_cellBClist.txt")

# 2) zUMIs 产生的bam文件
ZUMIS_BAM = os.path.join(OUT_DIR, f"{PROJECT_NAME}.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam")

# 3) bam_split 产生的文件
SPLIT_BAM_CELLTYPE1 = os.path.join(SPLIT_BAM_OUTPUT_DIR, f"{CELLTYPE1}_subset.bam")
SPLIT_BAM_CELLTYPE2 = os.path.join(SPLIT_BAM_OUTPUT_DIR, f"{CELLTYPE2}_subset.bam")

rule target:
    input:
        ZUMIS_BAM,
        SPLIT_BAM_CELLTYPE1,
        SPLIT_BAM_CELLTYPE2

rule extract_barcodes:
    input:
        barcode = ANNOTATION  
    output:
        CELLBC_CELLTYPE1,
        CELLBC_CELLTYPE2,
        CELLBC_CELLTYPE3
    params:
        script = os.path.join(algorithem_path, "SNP_calling", "extract_barcodes_and_split_bam.sh"),  
        output_folder = SPLIT_BAM_OUTPUT_DIR  
    log:
        "logs/extract_barcodelist.log"
    shell:
        """
        chmod +x {params.script}
        bash {params.script} {input.barcode} {params.output_folder}
        """

rule run_zUMIs:
    input:
        config_file
    params:
        script = os.path.join(ZUMIS_PATH, "zUMIs.sh")
    output:
        ZUMIS_BAM
    log:
        "logs/zUMIs.log"
    shell:
        "{params.script} -c -y {input}"

rule zUMIs_bam_index:
    input:
        ZUMIS_BAM
    output:
        ZUMIS_BAM + ".bai"
    shell:
        """
        samtools index {input}
        """

rule split_bam:
    input:
        bam = ZUMIS_BAM,
        BC_cell1 = CELLBC_CELLTYPE1,  
        BC_cell2 = CELLBC_CELLTYPE2  
    output:
        bam_cell1 = SPLIT_BAM_CELLTYPE1,
        bam_cell2 = SPLIT_BAM_CELLTYPE2
    params: 
        threads = THREADS
    shell:
        """
        samtools view -D BC:{input.BC_cell1} -@ {params.threads} -o {output.bam_cell1} {input.bam}
        samtools view -D BC:{input.BC_cell2} -@ {params.threads} -o {output.bam_cell2} {input.bam}
        """