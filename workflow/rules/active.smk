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
PROCESS_DIR = config["process_dir"]
MASKED_DIR = config["masked_dir"]

# ------------- #
# I/O           #
# ------------- #

# Masked genome for scanning
MASKED_GENOME = path.join(MASKED_DIR, "hg38", "mask", "hg38.custom-mask.fa")
MASKED_REGION = path.join(MASKED_DIR, "hg38", "mask", "hg38.masked_region.bed")

# Compiled matrix prob and scan
COMPILED_PROB = path.join("resources", "software", "PWMScan", "matrix_prob")
COMPILED_SCAN = path.join("resources", "software", "PWMScan", "matrix_scan")

# Premade list of JASPAR and UNIBND profiles
JASPAR_PROFILES_PATH = path.join("workflow", "misc", "profiles", "jaspar.txt")
UNIBIND_PROFILES_PATH = path.join("workflow", "misc", "profiles", "unibind.txt")

# Unibind TFBSs
UNIBIND_TFBS = path.join(
    "resources", "data", "unibind", "damo_hg38_TFBS_FLAT", "{profile}.bed"
)
UNIBIND_TFBS_MASKED = path.join(PROCESS_DIR, "mask", "tfbs", "{profile}.masked.bed")

# UniBind and JASPAR PFM download and unpacked output
JASPAR_MATRIX = path.join(
    INSTALL_DIR, "jaspar", "jaspar_pfms_redundant_FLAT", "{profile}.jaspar"
)
UNIBIND_FASTA = path.join(
    INSTALL_DIR, "unibind", "damo_hg38_FASTA_FLAT", "{profile}.fa"
)

# PFMS output
JASPAR_PFM = path.join(PROCESS_DIR, "pwm", "pfms", "jaspar", "{profile}", "counts.pfm")
UNIBIND_PFM = path.join(
    PROCESS_DIR, "pwm", "pfms", "unibind", "{profile}", "counts.pfm"
)
GENERAL_PFM = path.join(
    PROCESS_DIR, "pwm", "pfms", "{source}", "{profile}", "counts.pfm"
)

# PSSM output
PSSM = path.join(PROCESS_DIR, "pwm", "pssm", "{source}", "{profile}")

# Pvals and cutoff
PVALS = path.join(PROCESS_DIR, "pwm", "pssm", "{source}", "{profile}", "pvals.txt")
COEFF = path.join(PROCESS_DIR, "pwm", "pssm", "{source}", "{profile}", "coeff.txt")

# Summary of genome sites
UNIBIND_SITES = path.join(
    PROCESS_DIR,
    "sites",
    "hg38",
    "{source}",
    "{profile}",
    "sites.masked.summary.txt",
)

# Summary of jaspar sites
JASPAR_SITES = path.join(
    PROCESS_DIR,
    "sites",
    "hg38",
    "{source}",
    "{profile}",
    "sites.masked.summary.txt",
)

ACTIVITY_MAP = path.join(
    PROCESS_DIR,
    "activity",
    "{source}",
    "{profile}",
    "activity-map.tsv",
)

# Activity plot
ACTIVITY_PLT = path.join(
    PROCESS_DIR,
    "activity",
    "{source}",
    "{profile}",
    "activity-map.png",
)

# Biosample map
BIOSAMPLE_MAP = path.join(
    PROCESS_DIR,
    "report",
    "biosample_map.tsv",
)

# ------------- #
# Params        #
# ------------- #

# Thresholds
WINDOW_SIZE = config["window_size"]
PVAL_THRESH = config["pval_thresh"]
RVAL_THRESH = config["rval_thresh"]

# Types of data to download from UniBind
DOWNLOAD_TYPES = ["PWMS", "TFBS", "FASTA"]

# A little cheating - Read premade list of JASPAR and UNIBND profiles
with open(JASPAR_PROFILES_PATH, "r") as f:
    JASPAR_PROFILES = f.read().splitlines()
with open(UNIBIND_PROFILES_PATH, "r") as f:
    UNIBIND_PROFILES = f.read().splitlines()

# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        BIOSAMPLE_MAP,
        expand(ACTIVITY_PLT, source="unibind", profile=UNIBIND_PROFILES[:10]),


rule unibind_pfms:
    input:
        UNIBIND_FASTA,
    output:
        UNIBIND_PFM,
    params:
        profile=lambda wc: wc.profile,
    log:
        stdout="workflow/logs/unibind_pfms_{profile}.stdout",
        stderr="workflow/logs/unibind_pfms_{profile}.stderr",
    conda:
        "../envs/scan.yaml"
    threads: 1
    script:
        "../scripts/pfms/unibind.py"


rule jaspar_pfms:
    input:
        JASPAR_MATRIX,
    output:
        JASPAR_PFM,
    params:
        profile=lambda wc: wc.profile,
    log:
        stdout="workflow/logs/jaspar_pfms_{profile}.stdout",
        stderr="workflow/logs/jaspar_pfms_{profile}.stderr",
    conda:
        "../envs/scan.yaml"
    script:
        "../scripts/pfms/jaspar.py"


rule pssms:
    input:
        GENERAL_PFM,
    output:
        logo=path.join(PSSM, "logo.png"),
        pssm_log2=path.join(PSSM, "pssm.log2"),
        pssm_intLogOdds=path.join(PSSM, "pssm.intLogOdds"),
        pssm_stats=path.join(PSSM, "pssm.stats"),
    log:
        stdout="workflow/logs/pssms_{profile}_{source}.stdout",
        stderr="workflow/logs/pssms_{profile}_{source}.stderr",
    conda:
        "../envs/scan.yaml"
    threads: 1
    script:
        "../scripts/pssm/pssm.py"


rule calculate_probabilities:
    message:
        """
        Generates pval distribution for given PWM. Remove percentage signs and sort the input.
        """
    input:
        pssm=path.join(PSSM, "pssm.intLogOdds"),
        matrix_prob=COMPILED_PROB,
    output:
        PVALS,
    log:
        stdout="workflow/logs/calculate_probabilities_{profile}_{source}.stdout",
        stderr="workflow/logs/calculate_probabilities_{profile}_{source}.stderr",
    conda:
        "../envs/scan.yaml"
    shell:
        "{input.matrix_prob} {input.pssm} | sed 's/%//g' | sort -k1n  > {output}"


rule calculate_cutoff:
    message:
        """
        Relative threshold is percentage of best PWM.
        For 80%, format as integer 80. Makes life easier.
        Second awk grabs the head -n1, avoids pipefail
        """
    input:
        rules.calculate_probabilities.output,
    output:
        COEFF,
    params:
        pthresh=PVAL_THRESH,
        rthresh=RVAL_THRESH,
    log:
        stdout="workflow/logs/calculate_cutoff_{profile}_{source}.stdout",
        stderr="workflow/logs/calculate_cutoff_{profile}_{source}.stderr",
    conda:
        "../envs/scan.yaml"
    shell:
        """
        awk '{{if($2 < {params.pthresh} && $3>={params.rthresh}) print $1}}' {input} |
        awk 'FNR == 1' > {output}
        """


rule mask_unibind:
    """
    Filters UniBind sites out of custom masked regions used in scanning
    """
    input:
        unibind=UNIBIND_TFBS,
        exclude=MASKED_REGION,
    output:
        temp(UNIBIND_TFBS_MASKED),
    conda:
        "../envs/scan.yaml"
    log:
        stdout="workflow/logs/mask_unibind_{profile}.stdout",
        stderr="workflow/logs/mask_unibind_{profile}.stderr",
    shell:
        "bedtools intersect -a {input.unibind} -b {input.exclude} -v > {output}"


rule summarize_ix:
    """
    Intersect UniBind motifs with genome wide scanned motifs.
    - Intersect -c counts overlaps.
    - Vawk command makes overlaps binary (0/1), every so often more than one overlap.
    """
    input:
        pssm=rules.pssms.output.pssm_intLogOdds,
        coeff=rules.calculate_cutoff.output,
        hg38=MASKED_GENOME,
        matrix_scan=COMPILED_SCAN,
        unibind_sites=rules.mask_unibind.output,
    output:
        UNIBIND_SITES,
    conda:
        "../envs/scan.yaml"
    log:
        stdout="workflow/logs/summarize_ix_{profile}_{source}.stdout",
        stderr="workflow/logs/summarize_ix_{profile}_{source}.stderr",
    threads: 1
    shell:
        """
        {input.matrix_scan} -m {input.pssm} -c $(cat {input.coeff}) {input.hg38} |
        bedtools intersect -a stdin -b {input.unibind_sites} -c |
        vawk '{{if ($7>0) {{print $5, 1}} else {{print $5, 0}} }}' |
        sort -k1n |
        bedtools groupby -g 1 -c 2 -o freqdesc > {output}
        """


rule map_activity:
    """
    Creates activity map between UniBind and reference motifs.
    """
    input:
        sites=rules.summarize_ix.output,
        pvals=rules.calculate_probabilities.output,
    output:
        ACTIVITY_MAP,
    params:
        window=WINDOW_SIZE,
    conda:
        "../envs/scan.yaml"
    log:
        stdout="workflow/logs/map_activity_{source}_{profile}.stdout",
        stderr="workflow/logs/map_activity_{source}_{profile}.stderr",
    script:
        "../scripts/activity/map.py"


rule plot_activity:
    """
    Makes indvidual data plot based off activity map
    """
    input:
        activity=rules.map_activity.output,
    output:
        plt=ACTIVITY_PLT,
    params:
        profile=lambda wc: wc.profile,
    conda:
        "../envs/plot.yaml"
    log:
        stdout="workflow/logs/plot_activity_{source}_{profile}.stdout",
        stderr="workflow/logs/plot_activity_{source}_{profile}.stderr",
    script:
        "../scripts/activity/plot.R"


rule summarize_jaspar_sites:
    """
    Intersect UniBind motifs with genome wide scanned motifs.
    - Intersect -c counts overlaps.
    - Vawk command makes overlaps binary (0/1), every so often more than one overlap.
    """
    input:
        pssm=rules.pssms.output.pssm_intLogOdds,
        coeff=rules.calculate_cutoff.output,
        hg38=MASKED_GENOME,
        matrix_scan=COMPILED_SCAN,
    output:
        JASPAR_SITES,
    conda:
        "../envs/scan.yaml"
    log:
        stdout="workflow/logs/summarize_jaspar_sites_{profile}_{source}.stdout",
        stderr="workflow/logs/summarize_jaspar_sites_{profile}_{source}.stderr",
    threads: 1
    shell:
        """
        {input.matrix_scan} -m {input.pssm} -c $(cat {input.coeff}) {input.hg38} |
        sort -k5n |
        bedtools groupby -g 5 -c 1 -o count > {output}
        """


rule final_report:
    input:
        expand(ACTIVITY_PLT, source="unibind", profile=UNIBIND_PROFILES[:10]),
        expand(JASPAR_SITES, source="jaspar", profile=JASPAR_PROFILES[:10]),
    output:
        BIOSAMPLE_MAP,
    params:
        sites_dir=path.join("results", "unibind", "sites", "hg38"),
        pssm_dir=path.join("results", "unibind", "pwm", "pssm"),
        unibind_pwmdir=path.join("resources", "data", "unibind", "damo_hg38_PWMS_FLAT"),
        jaspar_profiles=path.join("workflow", "misc", "profiles", "jaspar.txt"),
        unibind_profiles=path.join("workflow", "misc", "profiles", "unibind.txt"),
        ambrosini_etal=path.join(
            "resources", "data", "ambrosini", "best4exp_all_jaspar.txt"
        ),
        activity_dir=path.join("results", "unibind", "activity"),
    log:
        stdout="workflow/logs/final_report.stdout",
        stderr="workflow/logs/final_report.stderr",
    conda:
        "../envs/scan.yaml"
    script:
        "../scripts/activity/report.py"
