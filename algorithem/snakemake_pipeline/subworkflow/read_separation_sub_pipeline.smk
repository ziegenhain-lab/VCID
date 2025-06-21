### read_separation_sub_pipeline.smk 子工作流 ----------------------------------------
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
CELLBCs_FOR_READSEPARATION = batch_config['readseparation_barcode_list']
Differetial_SNPs_from_other_source = config['Differetial_SNPs_from_other_source']


Data_results_path = batch_config['out_dir']
OUT_DIR = Data_results_path

# 3. Output files  ----------------------------------------------------------------
ZUMIS_OUTPUT_DIR = os.path.join(OUT_DIR, "zUMIs_output")

SPLIT_BAM_OUTPUT_DIR = os.path.join(OUT_DIR, "split_bam_output")

READ_SEPARATION_RESULTS_DIR = os.path.join(OUT_DIR, "read_separation_results")


# 4.Rule生成的文件的信息标注  -----------------------------------------------------
rule doublets_separation:
    input:
        yaml = config_file,
        vcf = (os.path.join(Data_results_path, "SNP_calling_output", f"{PROJECT_NAME}.processed.vcf")
               if not Differetial_SNPs_from_other_source else Differetial_SNPs_from_other_source),
        bam = os.path.join(OUT_DIR, f"{PROJECT_NAME}.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam"),
        barcodelist = CELLBCs_FOR_READSEPARATION,
        annotation = ANNOTATION
    params:
        subfolder = "read_separation_results", 
        script = os.path.join(algorithem_path, "DoubletSeparation_SNPbased_v1.0", "DoubletsSeparation_SNPbased_v0.6.0_separation.R")
    output:
        tmp_rds_path = os.path.join(READ_SEPARATION_RESULTS_DIR, f"{PROJECT_NAME}.tmp.rds"),
        readcount_matrix_rds_path = os.path.join(READ_SEPARATION_RESULTS_DIR, f"{PROJECT_NAME}.separated_counts.rds"),
        CELLTYPE1_readID = os.path.join(READ_SEPARATION_RESULTS_DIR, f"{PROJECT_NAME}.{CELLTYPE1}_readIDs.txt"),
        CELLTYPE2_readID = os.path.join(READ_SEPARATION_RESULTS_DIR, f"{PROJECT_NAME}.{CELLTYPE2}_readIDs.txt")
    shell:
        """
        Rscript {params.script} --yaml {input.yaml} --vcf {input.vcf} --bam {input.bam} \
        -c {input.barcodelist} -f {params.subfolder} -a {input.annotation}
        """

rule doublets_separation_plotting:
    input:
        rds = os.path.join(READ_SEPARATION_RESULTS_DIR, f"{PROJECT_NAME}.tmp.rds"),
        annotation = ANNOTATION
    params:
        output_dir = READ_SEPARATION_RESULTS_DIR,
        script = os.path.join(algorithem_path, "DoubletSeparation_SNPbased_v1.0", 
                        "DoubletsSeparation_SNPbased_v0.6.0_plotting.R")
    output:
        read_assignment_result = os.path.join(READ_SEPARATION_RESULTS_DIR, "final_read_assignment_result_sorted.txt")
    shell:
        """
        Rscript {params.script} --rds {input.rds} --annotation {input.annotation} \
        --outdir {params.output_dir} 
        """      

rule read_separation_for_bam:
    input:
        CELLTYPE1_readID = os.path.join(READ_SEPARATION_RESULTS_DIR, f"{PROJECT_NAME}.{CELLTYPE1}_readIDs.txt"),
        CELLTYPE2_readID = os.path.join(READ_SEPARATION_RESULTS_DIR, f"{PROJECT_NAME}.{CELLTYPE2}_readIDs.txt"),
        bam = os.path.join(OUT_DIR, f"{PROJECT_NAME}.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam")
    params:
        output_dir = READ_SEPARATION_RESULTS_DIR,
        script = os.path.join(algorithem_path, "DoubletSeparation_SNPbased_v1.0", "parallel_extract_bam.sh")
    output:
        celltype1_separated_bam_path = os.path.join(READ_SEPARATION_RESULTS_DIR, f"{PROJECT_NAME}.{CELLTYPE1}_readIDs_separated_bam.bam"),
        celltype2_separated_bam_path = os.path.join(READ_SEPARATION_RESULTS_DIR, f"{PROJECT_NAME}.{CELLTYPE2}_readIDs_separated_bam.bam")
    shell:
        """
        chmod +x {params.script}
        {params.script} {input.CELLTYPE1_readID} {input.bam} {params.output_dir} 50
        {params.script} {input.CELLTYPE2_readID} {input.bam} {params.output_dir} 50
        """      