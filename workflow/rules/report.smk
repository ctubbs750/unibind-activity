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


# ------------- #
# I/O           #
# ------------- #

# Input directories that need to be completed
PSSM_DIR = path.join(PROCESS_DIR, "pssm")
SCAN_DIR = path.join(PROCESS_DIR, "scan")

# Biosample map
BIOSAMPLE_MAP = path.join(PROCESS_DIR, "biosample_map.tsv")

# ------------- #
# Params        #
# ------------- #


# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
       BIOSAMPLE_MAP


rule make_mapping:
    input:
        PSSM_DIR,
        SCAN_DIR
    output:
        BIOSAMPLE_MAP,
    log:
        stdout="workflow/logs/make_mapping.stdout",
        stderr="workflow/logs/make_mapping.stderr",
    conda:
        "../envs/scan.yaml"
    script:
        "../scripts/report/mapping.py"