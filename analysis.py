"""
analysis.py

Reusable data-access and analysis functions used by both load_data.py's
downstream consumers and dashboard.py. All functions read from the SQLite
database produced by load_data.py.
"""

import sqlite3

import numpy as np
import pandas as pd
from scipy import stats

DB_FILE = "cell_count.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def get_connection(db_file=DB_FILE):
    return sqlite3.connect(db_file)


# ---------------------------------------------------------------------------
# Part 2: relative frequency summary table
# ---------------------------------------------------------------------------

def get_full_sample_table(conn=None):
    """Return one row per sample with subject/sample metadata attached."""
    own_conn = conn is None
    conn = conn or get_connection()
    query = """
        SELECT
            s.sample_id,
            s.subject_id,
            s.sample_type,
            s.time_from_treatment_start,
            sub.project,
            sub.condition,
            sub.age,
            sub.sex,
            sub.treatment,
            sub.response
        FROM samples s
        JOIN subjects sub ON s.subject_id = sub.subject_id
    """
    df = pd.read_sql_query(query, conn)
    if own_conn:
        conn.close()
    return df


def get_frequency_table(conn=None):
    """
    Part 2: for each sample, compute total cell count and the relative
    frequency (%) of each of the five populations.

    Returns a long-format DataFrame with columns:
        sample, total_count, population, count, percentage
    """
    own_conn = conn is None
    conn = conn or get_connection()

    counts = pd.read_sql_query(
        "SELECT sample_id AS sample, population, count FROM cell_counts", conn
    )
    if own_conn:
        conn.close()

    totals = counts.groupby("sample")["count"].sum().rename("total_count")
    counts = counts.merge(totals, on="sample")
    counts["percentage"] = (counts["count"] / counts["total_count"] * 100).round(3)

    counts = counts[["sample", "total_count", "population", "count", "percentage"]]
    counts = counts.sort_values(["sample", "population"]).reset_index(drop=True)
    return counts


# ---------------------------------------------------------------------------
# Part 3: responders vs non-responders statistical analysis
# ---------------------------------------------------------------------------

def get_miraclib_melanoma_pbmc_frequencies(conn=None):
    """
    Build the responder-vs-non-responder analysis dataset: relative
    frequencies for melanoma PBMC samples from subjects treated with
    miraclib, restricted to subjects with a recorded response.
    """
    own_conn = conn is None
    conn = conn or get_connection()

    freq = get_frequency_table(conn)
    meta = get_full_sample_table(conn)
    if own_conn:
        conn.close()

    merged = freq.merge(meta, left_on="sample", right_on="sample_id")

    mask = (
        (merged["condition"].str.lower() == "melanoma")
        & (merged["sample_type"].str.upper() == "PBMC")
        & (merged["treatment"].str.lower() == "miraclib")
        & (merged["response"].isin(["yes", "no"]))
    )
    return merged.loc[mask].reset_index(drop=True)


def get_subject_level_frequencies(conn=None):
    """
    Same cohort as get_miraclib_melanoma_pbmc_frequencies(), but collapsed
    to one row per (subject, population) by averaging across each subject's
    samples. Used for both the boxplot and the significance test so a
    subject with more timepoints in the data doesn't get more weight than
    one with fewer.
    """
    data = get_miraclib_melanoma_pbmc_frequencies(conn)
    return (
        data.groupby(["subject_id", "population", "response"])["percentage"]
        .mean()
        .reset_index()
    )


def compare_responders_vs_nonresponders(conn=None, alpha=0.05):
    """
    Part 3: for each population, compare relative-frequency distributions
    between responders and non-responders using the Mann-Whitney U test
    (non-parametric, appropriate for small/uneven clinical group sizes and
    doesn't assume normality). P-values are adjusted for multiple testing
    across the 5 populations using Benjamini-Hochberg FDR correction.

    Note on unit of analysis: subjects in this dataset contribute multiple
    samples each (one per timepoint). Testing at the sample level would
    treat those repeated samples as independent observations and inflate n,
    which isn't right since "responder" is a property of the patient, not
    the sample. So each subject's samples are first averaged per
    population, giving one value per subject before running the test.

    Returns a DataFrame with one row per population.
    """
    per_subject = get_subject_level_frequencies(conn)

    results = []
    for pop in POPULATIONS:
        sub = per_subject[per_subject["population"] == pop]
        resp = sub.loc[sub["response"] == "yes", "percentage"]
        nonresp = sub.loc[sub["response"] == "no", "percentage"]

        if len(resp) < 2 or len(nonresp) < 2:
            results.append({
                "population": pop,
                "n_responders": len(resp),
                "n_non_responders": len(nonresp),
                "responder_mean_pct": resp.mean() if len(resp) else np.nan,
                "non_responder_mean_pct": nonresp.mean() if len(nonresp) else np.nan,
                "u_statistic": np.nan,
                "p_value": np.nan,
            })
            continue

        u_stat, p_val = stats.mannwhitneyu(resp, nonresp, alternative="two-sided")
        results.append({
            "population": pop,
            "n_responders": len(resp),
            "n_non_responders": len(nonresp),
            "responder_mean_pct": round(resp.mean(), 3),
            "non_responder_mean_pct": round(nonresp.mean(), 3),
            "u_statistic": round(u_stat, 3),
            "p_value": p_val,
        })

    result_df = pd.DataFrame(results)

    # Benjamini-Hochberg FDR correction across the tested populations.
    valid = result_df["p_value"].notna()
    result_df["p_adj"] = np.nan
    if valid.sum() > 0:
        pvals = result_df.loc[valid, "p_value"].values
        order = np.argsort(pvals)
        ranked = pvals[order]
        m = len(ranked)
        adj = ranked * m / (np.arange(m) + 1)
        # enforce monotonicity
        adj = np.minimum.accumulate(adj[::-1])[::-1]
        adj = np.clip(adj, 0, 1)
        out = np.empty(m)
        out[order] = adj
        result_df.loc[valid, "p_adj"] = out

    result_df["significant"] = result_df["p_adj"] < alpha
    return result_df.sort_values("p_value").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Part 4: baseline miraclib melanoma PBMC subset analysis
# ---------------------------------------------------------------------------

def get_baseline_miraclib_melanoma_pbmc_samples(conn=None):
    """
    Part 4.1: melanoma PBMC samples at baseline (time_from_treatment_start
    == 0) from subjects treated with miraclib.
    """
    own_conn = conn is None
    conn = conn or get_connection()
    meta = get_full_sample_table(conn)
    if own_conn:
        conn.close()

    mask = (
        (meta["condition"].str.lower() == "melanoma")
        & (meta["sample_type"].str.upper() == "PBMC")
        & (meta["treatment"].str.lower() == "miraclib")
        & (meta["time_from_treatment_start"] == 0)
    )
    return meta.loc[mask].reset_index(drop=True)


def summarize_baseline_subset(conn=None):
    """
    Part 4.2: among the baseline subset, summarize:
      - sample counts per project
      - responder / non-responder counts
      - male / female counts
    Returns a dict of three DataFrames: by_project, by_response, by_sex.
    """
    subset = get_baseline_miraclib_melanoma_pbmc_samples(conn)

    by_project = (
        subset.groupby("project")["sample_id"]
        .count()
        .rename("n_samples")
        .reset_index()
        .sort_values("n_samples", ascending=False)
    )

    by_response = (
        subset.drop_duplicates("subject_id")
        .groupby("response")["subject_id"]
        .count()
        .rename("n_subjects")
        .reset_index()
    )

    by_sex = (
        subset.drop_duplicates("subject_id")
        .groupby("sex")["subject_id"]
        .count()
        .rename("n_subjects")
        .reset_index()
    )

    return {"by_project": by_project, "by_response": by_response, "by_sex": by_sex, "subset": subset}
