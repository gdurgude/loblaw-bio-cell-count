#!/usr/bin/env python3
"""Run the complete Loblaw Bio analysis pipeline."""

from pathlib import Path
import sqlite3
import pandas as pd
import plotly.express as px

import load_data
import analysis

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "output"
DB = ROOT / "cell_count.db"
CSV = ROOT / "cell-count.csv"


def save_plot(fig, filename):
    fig.write_html(OUTPUT / filename, include_plotlyjs="cdn")


def main():
    OUTPUT.mkdir(exist_ok=True)

    # Part 1: always rebuild the database from the supplied input.
    load_data.main()

    conn = sqlite3.connect(DB)
    try:
        # Part 2
        freq = analysis.get_frequency_table(conn)
        freq.to_csv(OUTPUT / "part2_frequency_table.csv", index=False)

        avg_comp = (
            freq.groupby("population", as_index=False)["percentage"]
            .mean()
            .rename(columns={"percentage": "mean_percentage"})
        )
        avg_comp.to_csv(OUTPUT / "part2_average_composition.csv", index=False)
        save_plot(
            px.bar(
                avg_comp,
                x="population",
                y="mean_percentage",
                labels={"mean_percentage": "Mean % of total cells"},
                title="Average cell-population composition",
            ),
            "part2_average_composition.html",
        )

        # Part 3
        stats = analysis.compare_responders_vs_nonresponders(conn)
        stats.to_csv(OUTPUT / "part3_statistical_results.csv", index=False)

        subject_freq = analysis.get_subject_level_frequencies(conn)
        subject_freq.to_csv(OUTPUT / "part3_subject_level_frequencies.csv", index=False)
        save_plot(
            px.box(
                subject_freq,
                x="population",
                y="percentage",
                color="response",
                points="all",
                category_orders={"response": ["yes", "no"]},
                labels={"percentage": "Relative frequency (%)"},
                title="Responders vs non-responders",
            ),
            "part3_responders_vs_nonresponders.html",
        )

        # Part 4
        summary = analysis.summarize_baseline_subset(conn)
        summary["subset"].to_csv(OUTPUT / "part4_baseline_subset.csv", index=False)
        summary["by_project"].to_csv(OUTPUT / "part4_by_project.csv", index=False)
        summary["by_response"].to_csv(OUTPUT / "part4_by_response.csv", index=False)
        summary["by_sex"].to_csv(OUTPUT / "part4_by_sex.csv", index=False)

        # Explicit requested answers.
        all_sex = analysis.count_subjects_by_sex(conn)
        baseline_sex = analysis.baseline_subject_counts_by_sex(conn)
        bcell = analysis.melanoma_male_responder_baseline_b_cell_average(conn)

        answers = [
            ["full_dataset_male_subjects", int(all_sex.loc[all_sex.sex.str.upper() == "M", "n_subjects"].sum()), ""],
            ["full_dataset_female_subjects", int(all_sex.loc[all_sex.sex.str.upper() == "F", "n_subjects"].sum()), ""],
            ["baseline_miraclib_melanoma_pbmc_male_subjects",
             int(baseline_sex.loc[baseline_sex.sex.str.upper() == "M", "n_subjects"].sum()), ""],
            ["baseline_miraclib_melanoma_pbmc_female_subjects",
             int(baseline_sex.loc[baseline_sex.sex.str.upper() == "F", "n_subjects"].sum()), ""],
            ["melanoma_male_responders_time0_average_b_cells",
             f"{float(bcell.iloc[0]['average_b_cells']):.2f}",
             f"n={int(bcell.iloc[0]['n_samples'])} baseline samples; all sample/treatment types"],
        ]
        pd.DataFrame(answers, columns=["metric", "value", "calculation_note"]).to_csv(
            OUTPUT / "requested_answers.csv", index=False
        )
    finally:
        conn.close()

    print(f"Pipeline complete. Outputs written to {OUTPUT}")


if __name__ == "__main__":
    main()
