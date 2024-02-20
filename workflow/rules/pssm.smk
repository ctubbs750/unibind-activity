from os import listdir, path
from pathlib import Path
from snakemake.utils import min_version

# Settings
min_version("7.32.4")


# ------------- #
# Config        #
# ------------- #

# I/O directories
INSTALL_DIR = config["install_dir"]
PROCESS_DIR = config["process_dir"]

# Download URLs
UNIBIND_URLS = {
    "PWMS": config["urls"]["pwms"],
    "TFBS": config["urls"]["tfbs"],
    "FASTA": config["urls"]["fasta"],
}

# JASPAR URL
JASPAR_PFMS_URL = config["urls"]["jaspar_pfms"]


# ------------- #
# I/O           #
# ------------- #

# Raw UniBind download and unpacked, flattened output
UNIBIND_DOWNLOAD = path.join(INSTALL_DIR, "damo_hg38_{type}.tar.gz")
UNIBIND_UNPACKED = path.join(INSTALL_DIR, "damo_hg38_{type}_FLAT")

# PFMS output
JASPAR_PFMS = path.join(PROCESS_DIR, "pwm","pfms", "jaspar")
DAMO_PFMS = path.join(PROCESS_DIR, "pwm", "pfms", "damo")

# PSSM output
JASPAR_PSSM = path.join(PROCESS_DIR, "motifs","pssm", "jaspar")
DAMO_PSSM = path.join(PROCESS_DIR, "motifs","pssm", "damo")

# Biosample maps
BIOSAMPLE_MAP = path.join(PROCESS_DIR, "biosample_map.tsv")

# ------------- #
# Params        #
# ------------- #


wildcard_constraints:
    tf_name="\w+",
    profile="MA\d{4}\.\d",
    biosample="[^/]+",


# Types of data to download from UniBind
DOWNLOAD_TYPES = ["PWMS", "TFBS", "FASTA"]

# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        # BIOSAMPLE_MAP,
        expand(UNIBIND_UNPACKED, type=DOWNLOAD_TYPES),
        JASPAR_PSSM,
        DAMO_PSSM,

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


rule make_damo_pfms:
    input:
        expand(rules.unpack_unibind.output, type=["FASTA"]),
    output:
        directory(DAMO_PFMS),
    log:
        stdout="workflow/logs/make_damo_pfms.stdout",
        stderr="workflow/logs/make_damo_pfms.stderr",
    conda:
        "../envs/scan.yaml"
    threads: 24
    script:
        "../scripts/pfms/damo.py"


rule make_jaspar_pfms:
    input:
        expand(rules.unpack_unibind.output, type=["FASTA"]),
    output:
        directory(JASPAR_PFMS),
    params:
        url=JASPAR_PFMS_URL,
    log:
        stdout="workflow/logs/make_jaspar_pfms.stdout",
        stderr="workflow/logs/make_jaspar_pfms.stderr",
    conda:
        "../envs/scan.yaml"
    script:
        "../scripts/pfms/jaspar.py"


rule make_damo_pssms:
    input:
        rules.make_damo_pfms.output,
    output:
        directory(DAMO_PSSM),
    log:
        stdout="workflow/logs/make_pssms.stdout",
        stderr="workflow/logs/make_pssms.stderr",
    conda:
        "../envs/scan.yaml"
    threads: 40
    script:
        "../scripts/pssm/pssm.py"


rule make_jaspar_pssms:
    input:
        rules.make_jaspar_pfms.output,
    output:
        directory(JASPAR_PSSM),
    log:
        stdout="workflow/logs/make_jaspar_pssms.stdout",
        stderr="workflow/logs/make_jaspar_pssms.stderr",
    conda:
        "../envs/scan.yaml"
    threads: 40
    script:
        "../scripts/pssm/pssm.py"