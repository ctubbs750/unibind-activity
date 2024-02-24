from pathlib import Path
from os import listdir, path
from snakemake.utils import min_version

# Settings
min_version("7.32.4")

# ------------- #
# Config        #
# ------------- #

# I/O directories
INSTALL_DIR = config["install_dir"]
GENOME_DIR = config["genome_dir"]
GENCODE_DIR = config["gencode_dir"]
MASKED_DIR = config["masked_dir"]
PWWMSCAN_DIR = config["pwmscan_dir"]

# UniBind URLs
UNIBIND_URLS = {
    "PWMS": config["urls"]["pwms"],
    "TFBS": config["urls"]["tfbs"],
    "FASTA": config["urls"]["fasta"],
}

# Jaspar PFM URL
JASPAR_PFMS_URL = config["urls"]["jaspar_pfms"]

# Ambrosini URL
AMBROSINI_URL = config["urls"]["ambrosini"]

# PWMScan
MATRIX_PROB = config["matrix_prob"]
MATRIX_SCAN = config["matrix_scan"]

# ------------- #
# I/O           #
# ------------- #

# Ambrosini AUROC download
AMBROSINI_DOWNLOAD = path.join(INSTALL_DIR, "unibind", "ambrosini", "best4exp_all_jaspar.txt")

# Genome file
GENOME_FA_GZ = path.join(GENOME_DIR, "hg38", "hg38.fa.gz")
GENOME_FA_UZ = path.join(GENOME_DIR, "hg38", "hg38.fa")

# Mask genome files
BLACKLIST = os.path.join(GENOME_DIR, "hg38", "hg38.blacklist.bed")
UCSC_GAPS = os.path.join(GENOME_DIR, "hg38", "hg38.gaps.bed")
EXON_SEQS = os.path.join(GENCODE_DIR, "hg38", "gencode.hg38.exons.protein_coding.bed")

# Masked genome files
MASKED_REGION = path.join(MASKED_DIR, "hg38", "mask", "hg38.masked_region.bed")
MASKED_GENOME = path.join(MASKED_DIR, "hg38", "mask", "hg38.custom-mask.fa")

# Raw UniBind download and unpacked, flattened output
UNIBIND_DOWNLOAD = path.join(INSTALL_DIR, "unibind", "damo_hg38_{type}.tar.gz")
UNIBIND_UNPACKED = path.join(INSTALL_DIR, "unibind", "damo_hg38_{type}_FLAT")

# JASPAR PFM download and unpacked output
JASPAR_DOWNLOAD = path.join(INSTALL_DIR, "jaspar", "jaspar_pfms_redundant.zip")
JASPAR_UNPACKED = path.join(INSTALL_DIR, "jaspar", "jaspar_pfms_redundant_FLAT")

# Compiled matrix prob and scan
COMPILED_PROB = path.join(PWWMSCAN_DIR, "matrix_prob")
COMPILED_SCAN = path.join(PWWMSCAN_DIR, "matrix_scan")

# ------------- #
# Params        #
# ------------- #

# Types of data to download from UniBind
DOWNLOAD_TYPES = ["PWMS", "TFBS", "FASTA"]


# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        AMBROSINI_DOWNLOAD,
        MASKED_GENOME,
        COMPILED_PROB,
        COMPILED_SCAN,
        JASPAR_UNPACKED,
        expand(UNIBIND_UNPACKED, type=DOWNLOAD_TYPES),


# ------------- #
# Data download #
# ------------- #


rule download_abrosini:
    message:
        "Downloads AUROC from Ambrosini et al. 2020"
    output:
        AMBROSINI_DOWNLOAD,
    params:
        url=AMBROSINI_URL,
    log:
        stdout="workflow/logs/download_abrosini.stdout",
        stderr="workflow/logs/download_abrosini.stderr",
    conda:
        "../envs/scan.yaml"
    shell:
        "curl -o {output} {params.url}"


rule download_unibind:
    message:
        "Downloads Unibind PWMs and TFBS from the UniBind database"
    output:
        UNIBIND_DOWNLOAD,
    params:
        url=lambda wc: UNIBIND_URLS[wc.type],
    log:
        stdout="workflow/logs/download_unibind_{type}.stdout",
        stderr="workflow/logs/download_unibind_{type}.stderr",
    conda:
        "../envs/scan.yaml"
    shell:
        "curl -o {output} {params.url}"


rule download_jaspar:
    message:
        "Downloads JASPAR redundant PFM set"
    output:
        JASPAR_DOWNLOAD,
    params:
        url=JASPAR_PFMS_URL,
    log:
        stdout="workflow/logs/download_jaspar.stdout",
        stderr="workflow/logs/download_jaspar.stderr",
    conda:
        "../envs/scan.yaml"
    shell:
        "curl -o {output} {params.url}"


rule unpack_jaspar:
    message:
        "Unpacks JASPAR PFM download"
    input:
        rules.download_jaspar.output,
    output:
        directory(JASPAR_UNPACKED),
    log:
        stdout="workflow/logs/unpack_jaspar.stdout",
        stderr="workflow/logs/unpack_jaspar.stderr",
    conda:
        "../envs/scan.yaml"
    shell:
        "unzip -q {input} -d {output}"


rule unpack_unibind:
    message:
        "Unpacks downloaded archives into flattened directories"
    input:
        rules.download_unibind.output,
    output:
        directory(UNIBIND_UNPACKED),
    log:
        stdout="workflow/logs/unpack_unibind_{type}.stdout",
        stderr="workflow/logs/unpack_unibind_{type}.stderr",
    conda:
        "../envs/scan.yaml"
    shell:
        "mkdir -p {output} && tar --strip-components 2 -xzf {input} -C {output}"


# ------------- #
# Genome setup  #
# ------------- #


rule decompress_genome:
    message:
        "Decompresses reference genome for scanning"
    input:
        GENOME_FA_GZ,
    output:
        GENOME_FA_UZ,
    log:
        stdout="workflow/logs/decompress_genome.stdout",
        stderr="workflow/logs/decompress_genome.stderr",
    conda:
        "../envs/scan.yaml"
    shell:
        "gunzip {input} -c > {output}"


rule mask_regions:
    message:
        "Makes exclude regions in BED format for genome masking."
    input:
        exons=EXON_SEQS,
        blacklist=BLACKLIST,
    output:
        MASKED_REGION,
    params:
        gaps=UCSC_GAPS,
    log:
        stdout="workflow/logs/mask_regions.stdout",
        stderr="workflow/logs/mask_regions.stderr",
    conda:
        "../envs/scan.yaml"
    shell:
        """
        # Concatenate gaps, exons, and blacklist files
        cat {params.gaps} {input.exons} {input.blacklist} |
        # Print the first three columns of each line
        vawk '{{print $1, $2, $3}}' | 
        # Remove duplicate lines
        vawk '!seen[$1, $2, $3]++' |
        # Sort the output by the first two columns
        sort -k 1,1 -k2,2n > {output}
        """


rule mask_genome:
    message:
        "Creates custom-masked genome for scanning against."
    input:
        regions=rules.mask_regions.output,
        genome=rules.decompress_genome.output,
    output:
        MASKED_GENOME,
    log:
        stdout="workflow/logs/mask_genome.stdout",
        stderr="workflow/logs/mask_genome.stderr",
    conda:
        "../envs/scan.yaml"
    shell:
        "bedtools maskfasta -fi {input.genome} -bed {input.regions} -fo {output}"


# ------------- #
# PWMScan setup #
# ------------- #


rule compile_pwmscan:
    message:
        "Compiles PWMScan. Note c99 flag. Installs into resources/software by default."
    input:
        prob=MATRIX_PROB,
        scan=MATRIX_SCAN,
    output:
        compiled_prob=COMPILED_PROB,
        compiled_scan=COMPILED_SCAN,
    log:
        stdout="workflow/logs/compile_pwmscan.stdout",
        stderr="workflow/logs/compile_pwmscan.stderr",
    conda:
        "../envs/scan.yaml"
    shell:
        """
        gcc -std=c99 -o {output.compiled_prob} {input.prob} &&
        gcc -std=c99 -o {output.compiled_scan} {input.scan}
        """
