import logomaker
from Bio import motifs
from pandas import DataFrame
from matplotlib import pyplot as plt
from Bio.motifs.jaspar import calculate_pseudocounts

# Snakemake parameters
IP_PFM = snakemake.input[0]  # type: ignore
OP_PSSM_LOGO = snakemake.output.logo  # type: ignore
OP_PSSM_LOG2 = snakemake.output.pssm_log2  # type: ignore
OP_PSSM_intLogOdds = snakemake.output.pssm_intLogOdds  # type: ignore
OP_PSSM_STATS = snakemake.output.pssm_stats  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def read_pfm(filepath: str):
    """Read PFM from from file and return transposed matrix for Biopython"""
    with open(filepath) as handle:
        m = motifs.read(handle, "pfm")
    return m


def save_logo(pssm_data: DataFrame, outpath: str) -> None:
    """Save logo plot"""

    pssm_data.columns = ["A", "C", "G", "T"]
    # Mask all values below 0
    pssm_data = pssm_data.mask(pssm_data < 0, 0)

    # # Get index in 1-base
    pssm_data.index = [i + 1 for i in pssm_data.index]

    # Create Logo object
    pwm_logo = logomaker.Logo(pssm_data, color_scheme="classic", ax=None)

    # Style using Logo methods
    pwm_logo.style_spines(visible=False)
    pwm_logo.style_spines(spines=["left", "bottom"], visible=True)
    pwm_logo.style_xticks(rotation=0, fmt="%d", anchor=0)

    # Style using Axes methods
    pwm_logo.ax.set_ylabel("Bits", labelpad=-1)
    pwm_logo.ax.xaxis.set_ticks_position("none")
    pwm_logo.ax.xaxis.set_tick_params(pad=-1)
    pwm_logo.ax.set_title("d", color="r")

    # # Save figure
    plt.savefig(outpath, dpi=100)


def save_pssm(pfm: str) -> None:
    """Main program"""
    # Read input PFM
    pfm = read_pfm(pfm)

    # Calculate pseudocounts
    pfm.pseudocounts = calculate_pseudocounts(pfm)

    # Caclulate PSSM
    pssm = list(map(list, zip(*[pfm.pssm[nt] for nt in "ACGT"])))  # type: ignore

    # Grab IC and GC content
    ic = pfm.pssm.mean()
    gc = pfm.counts.gc_content

    # Make intlogodds PSSM
    intlogodds = (DataFrame.from_dict(pssm) * 100).round().astype(int)

    # Format ic and gc as stats file
    stats_df = DataFrame.from_dict({"IC": [ic], "GC": [gc]})

    # Write outputs
    save_logo((DataFrame.from_dict(pssm)), OP_PSSM_LOGO)
    (DataFrame.from_dict(pssm)).to_csv(
        OP_PSSM_LOG2, sep="\t", header=False, index=False
    )
    intlogodds.to_csv(OP_PSSM_intLogOdds, sep="\t", header=False, index=False)
    stats_df.to_csv(OP_PSSM_STATS, sep="\t", header=True, index=False)


def main() -> None:
    """Main program"""
    # Call main function
    save_pssm(IP_PFM)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
