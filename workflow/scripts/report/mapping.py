""""""

import coreapi
from pathlib import Path
from pandas import DataFrame, read_csv
import numpy as np
from os import path
from scipy.stats import pearsonr


# Snakemake parameters
PSSM_DIR = snakemake.input[0]  # type: ignore
SCAN_DIR = snakemake.input[1]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def get_profiles(
    ip_dir: str,
) -> list:
    """d"""
    return [i.name for i in Path(ip_dir).glob("*") if i.is_dir()]


def extract_jaspar_info(data: dict) -> dict:
    """Extracts jaspar info from api pulldown"""
    return {
        "species": data.get("species", [{"name": "NaN"}])[0]["name"],
        "class": data.get("class", ["NaN"])[0],
        "family": data.get("family", ["NaN"])[0],
        "source": data.get("source", "NaN"),
        "length": len(data.get("pfm", {"A": []})["A"]),
    }


def extract_unibind_info(data: dict) -> dict:
    """Extracts info from api pulldown"""
    # Parse results from data query
    return {
        "biological_condition": data.get("biological_condition", "NaN"),
        "jaspar_id": data.get("jaspar_id", "NaN")[0],
        "jaspar_version": int(data.get("tfbs")[0]["DAMO"][0]["jaspar_version"]),  # type: ignore
        "score_threshold": float(data.get("tfbs")[0]["DAMO"][0]["score_threshold"]),  # type: ignore
        "total_tfbs": int(data.get("tfbs")[0]["DAMO"][0]["total_tfbs"]),  # type: ignore
        "total_peaks": data.get("total_peaks", "NaN"),
    }


def return_numsites_summary(filepath: str) -> int:
    """d"""
    data = read_csv(filepath, header=None, sep="\s+", names=["count", "score"])
    return data["count"].sum()


def read_pssm(filepath: str) -> DataFrame:
    """Reads the PSSM file"""
    return read_csv(filepath, header=None, sep="\t", names=["A", "C", "G", "T"])


def return_pssm_stats(filepath: str) -> tuple:
    """Read PSSM stats from a file"""
    stats = read_csv(filepath, sep="\t")
    ic = stats["IC"][0]
    gc = stats["GC"][0]
    return ic, gc


def main() -> None:
    """Main program"""
    # Initialize a client
    client = coreapi.Client()

    # Load the schema documents
    jaspar_schema = client.get("https://jaspar.elixir.no/api/v1/docs")
    unibind_schema = client.get("https://unibind.uio.no/api/v1/docs")

    # Setup actions
    jaspar_action = ["matrix", "read"]
    unibind_action = ["datasets", "read"]

    # Get list of all damo and JASPAR profiles
    damo_profiles = get_profiles(path.join(PSSM_DIR, "damo"))

    profile_check = {}
    biosample_info = []
    # Interact with the API endpoint
    for profile in damo_profiles:
        print(profile)
        # Get the profile name
        jaspar_id = ".".join(profile.split(".")[-3:-1])
        # Get the TFname, cell type, exp
        name = profile.split(".")[2]
        cell = profile.split(".")[0]
        # Get the experiment
        epxn = profile.split(".")[1]
        ###
        # Unibind API
        ###
        # Interact with the API endpoint
        biosample_id = ".".join(profile.split(".")[:3])
        params = {"tf_id": biosample_id}
        result = client.action(unibind_schema, unibind_action, params=params)
        # Extract the data
        data = extract_unibind_info(result)
        data.update({"experiment": epxn, "biosample": cell, "name": name})
        ###
        # JASPAR API
        ###
        # Check if processed jaspar id
        if jaspar_id not in profile_check.keys():
            # Jaspar params
            params = {"matrix_id": jaspar_id}
            # Extract the data, handle failed query
            try:
                result = client.action(jaspar_schema, jaspar_action, params=params)
                jaspar = extract_jaspar_info(result)
            except:
                jaspar = {
                    "species": "NaN",
                    "class": "NaN",
                    "family": "NaN",
                    "source": "NaN",
                    "length": "NaN",
                }
            # Update profile check and data
            profile_check[profile] = jaspar
            data.update(jaspar)
            ###
            # Num sites from scanning -canonical
            ###
            # Paths
            scan_path = f"{path.join(SCAN_DIR, 'hg38', 'sites', 'jaspar', jaspar_id)}/sites.masked.summary.txt"
            # Get the number of sites
            try:
                data.update({"nsites_scan_canon": return_numsites_summary(scan_path)})
            except:
                data.update({"nsites_scan_canon": "NaN"})
        else:
            data.update(profile_check[profile])
        ###
        # Num sites from scanning
        ###

        # Paths
        scan_path = f"{path.join(SCAN_DIR, 'hg38', 'sites', 'damo', profile)}/sites.masked.summary.txt"
        # Get the number of sites
        data.update({"nsites_scan_damo": return_numsites_summary(scan_path)})

        ###
        # PCC
        ###

        # Get pearson correlation for canonical and original damo
        # Canonical
        try:
            pssm_path = f"{path.join(PSSM_DIR, 'jaspar', jaspar_id)}/pssm.intLogOdds"
            pssm = read_pssm(pssm_path)
            canon_vals_flat = pssm.values.ravel()

            # Damo
            pssm_path = f"{path.join(PSSM_DIR, 'damo', profile)}/pssm.intLogOdds"
            pssm = read_pssm(pssm_path)
            damo_vals_flat = pssm.values.ravel()

            # Original damo - ene
            pssm_path = path.join(
                "resources", "data", "unibind", "damo_hg38_PWMS_FLAT", f"{profile}.pwm"
            )
            pssm = read_csv(pssm_path, header=None, sep="\t").T
            damo_ene_flat = pssm.values.ravel()

            # Get PCCs
            data.update({"pcc_canon": pearsonr(canon_vals_flat, damo_vals_flat)[0]})
            data.update({"pcc_energy": pearsonr(damo_vals_flat, damo_ene_flat)[0]})
        except:
            data.update({"pcc_canon": "NaN"})
            data.update({"pcc_energy": "NaN"})

        ###
        # GC and IC
        ###
        stats_path = f"{path.join(PSSM_DIR, 'damo', profile)}/pssm.stats"
        gc, ic = return_pssm_stats(stats_path)

        # Update data
        data.update({"gc": gc, "ic": ic})

        # Store results
        biosample_info.append(data)

    # Convert to dataframe
    biosample_info = DataFrame(biosample_info)

    # Clean up biological condition - replace empty lists with NaN
    biosample_info["biological_condition"] = biosample_info[
        "biological_condition"
    ].apply(lambda x: np.nan if x == [] else x)

    # strip brackets from lists in biological condition
    biosample_info["biological_condition"] = biosample_info[
        "biological_condition"
    ].apply(lambda x: x[0] if type(x) == list else x)

    # Save to file
    biosample_info.to_csv(OUTPUT, index=False, sep="\t")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
