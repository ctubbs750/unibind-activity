""""""

import coreapi
import numpy as np
from os import path
from scipy.stats import pearsonr
from pandas import DataFrame, read_csv


# Snakemake parameters
SITES_DIR = snakemake.params["sites_dir"]  # type: ignore
PSSM_DIR = snakemake.params["pssm_dir"]  # type: ignore
UNIBIND_PWMDIR = snakemake.params["unibind_pwmdir"]  # type: ignore
JASPAR_PROFILES = snakemake.params["jaspar_profiles"]  # type: ignore
UNIBIND_PROFILES = snakemake.params["unibind_profiles"]  # type: ignore
AMBROSINI_ETAL = snakemake.params["ambrosini_etal"]  # type: ignore
ACTIVITY_DIR = snakemake.params["activity_dir"]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def get_profiles(filepath: str) -> DataFrame:
    """Reads in profiles"""
    return read_csv(filepath, header=None, names=["profile"], sep="\t", engine="c")[
        "profile"
    ]


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


def nscan_canon(filepath: str) -> int:
    """d"""
    # Readin
    data = read_csv(filepath, header=None, sep="\s+", names=["score", "count"])
    # Return sum of all sites
    return data["count"].sum()


def nscan_damo(filepath: str) -> int:
    """d"""
    # Readin
    data = read_csv(filepath, header=None, sep="\s+", names=["score", "summary"])
    # Test if both bound and unbound are present
    test = data["summary"].str.split(":").apply(lambda x: len(x)).max()

    if test > 2:
        # Get count summary info
        data[[0, 1, 2, 3]] = (
            data["summary"]
            .str.replace(",", ":")
            .str.split(":", expand=True)
            .fillna(0)
            .astype(int)
        )

        # Return sum of bound and unbound
        return data[1].sum() + data[2].sum()
    else:
        # Get count summary info
        data[[0, 1]] = (
            data["summary"]
            .str.replace(",", ":")
            .str.split(":", expand=True)
            .fillna(0)
            .astype(int)
        )
        # Return sum
        return data[1].sum()


def read_pssm(filepath: str) -> DataFrame:
    """Reads the PSSM file"""
    return read_csv(filepath, header=None, sep="\t", names=["A", "C", "G", "T"])


def return_pssm_stats(filepath: str) -> tuple:
    """Read PSSM stats from a file"""
    stats = read_csv(filepath, sep="\t")
    ic = stats["IC"][0]
    gc = stats["GC"][0]
    return ic, gc


def prepare_ambrosini(filepath: str) -> DataFrame:
    """Reads ambrosini et al. 2020 and makes dbase"""
    # Read in
    ambrosini = read_csv(
        filepath,
        sep="\t",
        header=None,
        names=["expn", "matrix", "tf_expn", "tf_matrix", "pv"],
    )
    # Get info
    ambrosini[["expn_id", "tf_name", "cell_tag1", "cell_tag2"]] = ambrosini[
        "expn"
    ].str.split(".", expand=True)
    # Reutrn dict
    return dict(list(zip(ambrosini["expn"], ambrosini["pv"])))


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
    unibind_profiles = get_profiles(UNIBIND_PROFILES)

    # Load ambrosini et al. 2020 database
    ambrosini_dict = prepare_ambrosini(AMBROSINI_ETAL)

    #
    profile_check = {}
    biosample_info = []
    # Interact with the API endpoint
    for profile in unibind_profiles:
        print(profile)
        # Get the profile name
        jaspar_id = ".".join(profile.split(".")[-3:-1])

        # Get the TFname, cell type, exp
        expn = profile.split(".")[0]
        cell = profile.split(".")[1]
        name = profile.split(".")[2]
        algo = profile.split(".")[5]

        ###
        # Unibind API
        ###

        # Make tf_id, Interact with the API endpoint
        tf_id = ".".join(profile.split(".")[:3])
        params = {"tf_id": tf_id}
        result = client.action(unibind_schema, unibind_action, params=params)

        # Extract the data
        data = extract_unibind_info(result)
        data.update({"expn": expn, "biosample": cell, "name": name, "algo": algo})

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
            scan_path = (
                f"{path.join(SITES_DIR, 'jaspar', jaspar_id)}/sites.masked.summary.txt"
            )

            # Get the number of sites
            try:
                data.update({"nscan_canon": nscan_canon(scan_path)})
            except:
                data.update({"nscan_canon": "NaN"})
        else:
            data.update(profile_check[profile])

        ###
        # Num sites from biosample scanning
        ###

        # Paths
        scan_path = (
            f"{path.join(SITES_DIR, 'unibind', profile)}/sites.masked.summary.txt"
        )

        # Get the number of sites
        data.update({"nsites_scan_damo": nscan_damo(scan_path)})

        ###
        # PCC - Canon and damo
        ###

        # Get pearson correlation for canonical and damo
        try:
            # Setup paths
            pssm_canon_path = (
                f"{path.join(PSSM_DIR, 'jaspar', jaspar_id)}/pssm.intLogOdds"
            )
            pssm_damo_log_path = f"{path.join(PSSM_DIR, 'unibind', profile)}/pssm.log2"
            pssm_damo_ene_path = path.join(UNIBIND_PWMDIR, f"{profile}.pwm")

            # Read in
            pssm_canon_data = read_pssm(pssm_canon_path)
            pssm_damo_log_data = read_pssm(pssm_damo_log_path)
            pssm_damo_ene_data = read_csv(pssm_damo_ene_path, header=None, sep="\t").T

            # Flatten
            canon_vals_flat = pssm_canon_data.values.ravel()
            damo_log_vals_flat = pssm_damo_log_data.values.ravel()
            damo_ene_vals_flat = pssm_damo_ene_data.values.ravel()

            # Get PCCs
            data.update({"pcc_canon": pearsonr(canon_vals_flat, damo_log_vals_flat)[0]})
            data.update(
                {"pcc_energy": pearsonr(damo_log_vals_flat, damo_ene_vals_flat)[0]}
            )
        except:
            data.update({"pcc_canon": "NaN"})
            data.update({"pcc_energy": "NaN"})

        ###
        # GC and IC
        ###

        stats_path = f"{path.join(PSSM_DIR, 'unibind', profile)}/pssm.stats"
        gc, ic = return_pssm_stats(stats_path)

        # Convert GC to proportion
        gc = gc / 100

        # Update data
        data.update({"gc": gc, "ic": ic})

        ###
        # Ambrosini et al. 2020 ROC AUC
        ###

        # Make key
        ambro_key = f"{expn}.{name}.{cell.split('_')[0]}"

        # Get AUC
        data.update({"auroc_ambro": ambrosini_dict.get(ambro_key, "NaN")})

        ###
        # Activity AUC
        ###

        # Get activity auc
        activity_path = path.join(
            ACTIVITY_DIR, "unibind", f"{profile}/activity-map.mse"
        )
        activity_mse = read_csv(activity_path, header=None, names=["mse"])["mse"].item()

        # Update
        data.update({"activity_mse": activity_mse})

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
