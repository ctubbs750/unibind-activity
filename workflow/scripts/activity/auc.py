"""d"""

from typing import List, Optional
from sklearn import metrics
from pandas import DataFrame, read_csv, merge
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split


# Snakemake
UNIBIND_IX = snakemake.input[0]  # type: ignore
PROFILE_PVALS = snakemake.input[1]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore
PROFILE = snakemake.params[0]  # type: ignore
DATASET = snakemake.params[1]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


# def read_unibind_ix(
#     filepath: str, fields: list = [4, 6], labels: list = ["score", "unibind"]
# ) -> DataFrame:
#     """Reads genome sites ix with unibind sites"""
#     # Return data
#     return read_csv(
#         filepath,
#         sep="\t",
#         header=None,
#         usecols=fields,
#         names=labels,
#         dtype=dict(zip(labels, [int, int])),
#     )


def read_unibind_ix(
    filepath: str, fields: List[int] = None, labels: List[str] = None
) -> DataFrame:
    """
    Reads genome sites ix with unibind sites

    Parameters:
    filepath (str): Path to the file to read
    fields (List[int]): List of column indices to use. Defaults to [4, 6].
    labels (List[str]): List of column names to use. Defaults to ["score", "unibind"].

    Returns:
    DataFrame: The read data
    """
    if fields is None:
        fields = [4, 6]
    if labels is None:
        labels = ["score", "unibind"]

    return read_csv(
        filepath,
        sep="\t",
        header=None,
        usecols=fields,
        names=labels,
        dtype=dict(zip(labels, [int, int])),
    )


# def read_profile_pvals(filepath: str, labels: list = ["score", "pval", "perc"]) -> dict:
#     """Reads pval data for provile"""
#     # Return data
#     return read_csv(
#         filepath,
#         delim_whitespace=True,
#         header=None,
#         names=labels,
#         dtype=dict(zip(labels, [int, float, float])),
#     )


def read_profile_pvals(filepath: str, labels: List[str] = None) -> DataFrame:
    """
    Reads pval data for profile

    Parameters:
    filepath (str): Path to the file to read
    labels (List[str]): List of column names to use. Defaults to ["score", "pval", "perc"].

    Returns:
    DataFrame: The read data
    """
    if labels is None:
        labels = ["score", "pval", "perc"]

    return read_csv(
        filepath,
        delim_whitespace=True,
        header=None,
        names=labels,
        dtype=dict(zip(labels, [int, float, float])),
    )


# def calculate_auc(
#     data: DataFrame, outcome: str = "unibind", covars: list = ["score"]
# ) -> float:
#     """Calculates AUC"""

#     # Setup outcome and covars
#     y = data[outcome]
#     X = data[covars]

#     # Split into training and testing
#     X_train, X_test, y_train, y_test = train_test_split(
#         X, y, test_size=0.3, random_state=0
#     )

#     # Fit logistic regression model
#     model = LogisticRegression(n_jobs=12, solver="saga")
#     print("fitting model...")
#     model.fit(X_train, y_train)

#     # Make predctions
#     print("making predictions...")
#     y_pred = model.predict_proba(X_test)[:, 1]


#     # Calculate and return auc
#     return round(metrics.roc_auc_score(y_test, y_pred), 4)
def calculate_auc(
    data: DataFrame,
    outcome: str = "unibind",
    covars: List[str] = None,
    random_state: int = 0,
) -> float:
    """
    Calculates AUC

    Parameters:
    data (DataFrame): The data to use
    outcome (str): The outcome variable. Defaults to "unibind".
    covars (List[str]): The covariate(s). Defaults to ["score"].
    random_state (int): The random state for train-test split. Defaults to 0.

    Returns:
    float: The calculated AUC
    """
    if covars is None:
        covars = ["score"]

    # Setup outcome and covars
    y = data[outcome]
    X = data[covars]

    # Split into training and testing
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=random_state
    )

    # Fit logistic regression model
    model = LogisticRegression(n_jobs=12, solver="saga")
    model.fit(X_train, y_train)

    # Make predictions
    y_pred = model.predict_proba(X_test)[:, 1]

    # Calculate and return auc
    return round(metrics.roc_auc_score(y_test, y_pred), 4)


# def main():
#     """Main program"""
#     # Read intersection and pvals data
#     unibind_ix = read_unibind_ix(UNIBIND_IX)
#     profile_pvals = read_profile_pvals(PROFILE_PVALS)

#     # Merge pvals with ix data to convert PWM scores to percentage
#     unibind_ix = merge(unibind_ix, profile_pvals, on="score", how="left")

#     # Calculate AUC - both raw PWM and perc
#     try:
#         auc_pwm = calculate_auc(data=unibind_ix, covars=["score"])
#         auc_perc = calculate_auc(data=unibind_ix, covars=["perc"])
#     except:
#         auc_pwm = "NA"
#         auc_perc = "NA"

#     # DataFrame
#     results = DataFrame(
#         {
#             "auc_pwm": [auc_pwm],
#             "auc_perc": [auc_perc],
#             "profile": [PROFILE],
#             "dataset": [DATASET],
#         }
#     )
#     # Write out
#     results.to_csv(OUTPUT, index=False, sep="\t", header=None)


def main(
    unibind_ix_file: str, profile_pvals_file: str, profile: str, dataset: str
) -> Optional[DataFrame]:
    """
    Main program

    Parameters:
    unibind_ix_file (str): Path to the unibind_ix file
    profile_pvals_file (str): Path to the profile_pvals file
    profile (str): The profile
    dataset (str): The dataset

    Returns:
    DataFrame: The results, or None if an error occurred
    """
    # Read intersection and pvals data
    unibind_ix = read_unibind_ix(unibind_ix_file)
    profile_pvals = read_profile_pvals(profile_pvals_file)

    # Merge pvals with ix data to convert PWM scores to percentage
    unibind_ix = merge(unibind_ix, profile_pvals, on="score", how="left")

    # Calculate AUC - both raw PWM and perc
    try:
        auc_pwm = calculate_auc(data=unibind_ix, covars=["score"])
        auc_perc = calculate_auc(data=unibind_ix, covars=["perc"])
    except:
        auc_pwm = "NA"
        auc_perc = "NA"

    # DataFrame
    results = DataFrame(
        {
            "auc_pwm": [auc_pwm],
            "auc_perc": [auc_perc],
            "profile": [profile],
            "dataset": [dataset],
        }
    )

    return results


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    results = main(UNIBIND_IX, PROFILE_PVALS, PROFILE, DATASET)
    if results is not None:
        results.to_csv(OUTPUT, index=False, sep="\t", header=None)
