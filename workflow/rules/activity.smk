from os import listdir, path
from pathlib import Path
from snakemake.utils import min_version

# Settings
min_version("7.32.4")


# ------------- #
# Config        #
# ------------- #

# FORMAT = config["format"]
ASSEMBLY = config["assembly"]
OUTPUT_DIR = config["output_dir"]
GENOME_DIR = config["genome_dir"]
GENCODE_DIR = config["gencode_dir"]
MATRIX_PROB = config["matrix_prob"]
MATRIX_SCAN = config["matrix_scan"]
CHROMOSOMES = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]

# NEW
PSSM_DIR = "results/unibind/pwm/pssm"

# ------------- #
# I/O           #
# ------------- #

# Unibind TFBSs
UNIBIND_SITES = path.join(
    "resources", "data", "unibind", "damo_hg38_TFBS_FLAT", "{profile}.bed"
)
UNIBIND_SITES_MASKED = path.join(
    "results", "unibind", "sites", "masked", "{profile}.masked.bed"
)

# Input PSSM structure
PSSM = path.join(PSSM_DIR, "{source}", "{profile}", "pssm.intLogOdds")

# Genome files
GENOME_FA = os.path.join(GENOME_DIR, ASSEMBLY, f"{ASSEMBLY}.fa.gz")
BLACKLIST = os.path.join(GENOME_DIR, ASSEMBLY, f"{ASSEMBLY}.blacklist.bed")
UCSC_GAPS = os.path.join(GENOME_DIR, ASSEMBLY, f"{ASSEMBLY}.gaps.bed")

# Gencode protein-coding exons
EXONS = os.path.join(GENCODE_DIR, ASSEMBLY, "gencode.hg38.exons.protein_coding.bed")

# Masked genome files
MASKED_SITES = path.join(OUTPUT_DIR, ASSEMBLY, "mask", f"{ASSEMBLY}.masked_sites.bed")
MASKED_GENOME = path.join(OUTPUT_DIR, ASSEMBLY, "mask", f"{ASSEMBLY}.custom-mask.fa")

# Compiled matrix prob and scan
COMPILED_PROB = os.path.join("resources", "software", "PWMScan", "matrix_prob")
COMPILED_SCAN = os.path.join("resources", "software", "PWMScan", "matrix_scan")

# Pvals and cutoff
PVALS = path.join(OUTPUT_DIR, "pssm", "{source}", "{profile}", "pvals.txt")
COEFF = path.join(OUTPUT_DIR, "pssm", "{source}", "{profile}", "coeff.txt")

# Summary of genome sites
SITES = path.join(
    OUTPUT_DIR,
    ASSEMBLY,
    "sites",
    "{source}",
    "{profile}",
    "sites.masked.summary.txt",
)
SITES_FOR_AUC = path.join(
    OUTPUT_DIR,
    ASSEMBLY,
    "sites",
    "{source}",
    "{profile}",
    "sites.for_auc.txt.gz",
)

SITES_FOR_AUC_PVALS = path.join(
    OUTPUT_DIR,
    ASSEMBLY,
    "sites",
    "{source}",
    "{profile}",
    "sites.for_auc.pvals.txt.gz",
)

ACTIVITY_MAP = os.path.join(
    "results",
    "unibind",
    "activity",
    "{source}",
    "{profile}",
    "activity-map.tsv",
)

# Activity plot
ACTIVITY_PLT = os.path.join(
    "results",
    "unibind",
    "activity",
    "{source}",
    "{profile}",
    "activity.map.png",
)

# ------------- #
# Params        #
# ------------- #

JASPAR_PROFILE_NAMES = [
    i.name for i in Path(path.join(PSSM_DIR, "jaspar")).iterdir() if i.is_dir()
]
DAMO_BIOSAMPLE_NAMES = [
    i.name for i in Path(path.join(PSSM_DIR, "damo")).iterdir() if i.is_dir()
]

# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        # expand(SITES, source="jaspar", profile=JASPAR_PROFILE_NAMES),
        # expand(SITES, source="damo", profile=DAMO_BIOSAMPLE_NAMES),
        # expand(SITES_FOR_AUC_PVALS, source="damo", profile=DAMO_BIOSAMPLE_NAMES),
        #expand(ACTIVITY_MAP, source="damo", profile=DAMO_BIOSAMPLE_NAMES)[:10],
        expand(ACTIVITY_PLT, source="damo", profile=DAMO_BIOSAMPLE_NAMES)[:10],


rule decompress_genome:
    message:
        "Necessary to scan genome."
    input:
        GENOME_FA,
    output:
        temp(f"{path.splitext(GENOME_FA)[0]}.fa"),
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
        exons=EXONS,
        blacklist=BLACKLIST,
    output:
        temp(MASKED_SITES),
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


rule compile_pwmscan:
    message:
        "Compiles PWMScan. Note c99 flag. Installs into resources/software by default."
    input:
        prob="../scripts/scan/matrix_prob.c",
        scan="../scripts/scan/matrix_scan.c",
    output:
        compiled_prob=COMPILED_PROB,
        compiled_scan=COMPILED_SCAN,
    log:
        stdout="workflow/logs/compile_pwmscan.stdout",
        stderr="workflow/logs/compile_pwmscan.stderr",
    conda:
        "../envs/scan.yaml"
    shell:
        "gcc -std=c99 -o {output.compiled_prob} {input.prob} && gcc -std=c99 -o {output.compiled_scan} {input.scan}"


rule calculate_probabilities:
    message:
        """
        Generates pval distribution for given PWM.
        # Remove percentage signs and sort the input
        """
    input:
        pssm=PSSM,
        matrix_prob=rules.compile_pwmscan.output.compiled_prob,
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
        pthresh=0.05,
        rthresh=90,
    log:
        stdout="workflow/logs/calculate_cutoff_{profile}_{source}.stdout",
        stderr="workflow/logs/calculate_cutoff_{profile}_{source}.stderr",
    conda:
        "../envs/scan.yaml"
    shell:
        "awk '{{if($2 < {params.pthresh} && $3>={params.rthresh}) print $1}}' {input} | awk 'FNR == 1' > {output}"


rule sites_summary:
    message:
        "Scans reference genome for matches against input motif. Output compressed."
    input:
        pssm=PSSM,
        coeff=rules.calculate_cutoff.output,
        hg38=rules.mask_genome.output,
        matrix_scan=rules.compile_pwmscan.output.compiled_scan,
    output:
        SITES,
    log:
        stdout="workflow/logs/sites_summary_{profile}_{source}.stdout",
        stderr="workflow/logs/sites_summary_{profile}_{source}.stderr",
    conda:
        "../envs/scan.yaml"
    shell:
        "{input.matrix_scan} -m {input.pssm} -c $(cat {input.coeff}) {input.hg38} | cut -f5 | sort | uniq -c > {output}"


rule mask_unibind:
    """
    Filters UniBind sites out of custom masked regions used in scanning
    """
    input:
        unibind=UNIBIND_SITES,
        exclude=MASKED_SITES,
    output:
        temp(UNIBIND_SITES_MASKED),
    conda:
        "../envs/scan.yaml"
    log:
        stdout="workflow/logs/mask_unibind_{profile}.stdout",
        stderr="workflow/logs/mask_unibind_{profile}.stderr",
    shell:
        "bedtools intersect -a {input.unibind} -b {input.exclude} -v > {output}"


rule intersect_motifs:
    """
    Intersect UniBind motifs with genome wide scanned motifs.
    - Intersect -c counts overlaps.
    - Vawk command makes overlaps binary (0/1), every so often more than one overlap.
    """
    input:
        pssm=PSSM,
        coeff=rules.calculate_cutoff.output,
        hg38=rules.mask_genome.output,
        matrix_scan=rules.compile_pwmscan.output.compiled_scan,
        unibind_sites=rules.mask_unibind.output,
    output:
        temp(SITES_FOR_AUC),
    conda:
        "../envs/scan.yaml"
    log:
        stdout="workflow/logs/intersect_motifs_{profile}_{source}.stdout",
        stderr="workflow/logs/intersect_motifs_{profile}_{source}.stderr",
    threads: 1
    shell:
        """
        {input.matrix_scan} -m {input.pssm} -c $(cat {input.coeff}) {input.hg38} |
        bedtools intersect -a stdin -b {input.unibind_sites} -c |
        vawk '{{if ($7>0) {{print $5, 1}} else {{print $5, 0}} }}' | gzip > {output}
        """


rule update_sites_wpvals:
    """
    Adds pvals to sites for AUC.
    """
    input:
        sites=rules.intersect_motifs.output,
        pvals=rules.calculate_probabilities.output,
    output:
        SITES_FOR_AUC_PVALS,
    conda:
        "../envs/scan.yaml"
    log:
        stdout="workflow/logs/update_sites_wpvals_{profile}_{source}.stdout",
        stderr="workflow/logs/update_sites_wpvals_{profile}_{source}.stderr",
    threads: 1
    script:
        "../scripts/scan/pvals.py"


rule map_activity:
    """
    Creates activity map between UniBind and reference motifs.
    """
    input:
        SITES_FOR_AUC_PVALS,
    output:
        ACTIVITY_MAP,
    params:
        # window=10,
        # threshold=0.10,
        window=10,
        threshold=10,
    conda:
        "../envs/scan.yaml"
    log:
        stdout="workflow/logs/map_activity_{source}_{profile}.stdout",
        stderr="workflow/logs/map_activity_{source}_{profile}.stderr",
    script:
        "../scripts/activity/activity.py"

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
