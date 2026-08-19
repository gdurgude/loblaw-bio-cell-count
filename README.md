# Loblaw Bio-Cell Count Analysis

This repository contains a reproducible Python/SQLite pipeline for the immune-cell count assignment. The pipeline loads the supplied CSV into a normalized relational database, produces the requested analyses and output files, and serves an interactive Streamlit dashboard.

## Quick start

The project is designed for GitHub Codespaces.

```bash
make setup
make pipeline
make dashboard
```

- `make setup` installs all Python dependencies.
- `make pipeline` rebuilds the SQLite database and runs Parts 1–4 from start to finish.
- `make dashboard` runs the Streamlit dashboard on port `8501`. In Codespaces, open the forwarded port 8501. Locally: http://localhost:8501

No manual database initialization or Python script execution is required.

## Repository structure

```text
.
├── cell-count.csv                 # supplied input data
├── load_data.py                   # database schema + CSV loader
├── analysis.py                    # reusable analysis/query functions
├── pipeline.py                    # complete end-to-end pipeline
├── dashboard.py                   # Streamlit dashboard
├── requirements.txt               # Python dependencies
├── Makefile                       # required grader entry points
├── output/                        # generated tables and interactive plots
└── README.md
```

## Database schema

The database is SQLite and is rebuilt by `make pipeline`.

The source CSV contains information at three different grains, so the database separates those grains:

```text
subjects
---------
subject_id       PK
project
condition
age
sex
treatment
response

        1
        |
        | many
        v
samples
-------
sample_id                    PK
subject_id                   FK -> subjects.subject_id
sample_type
time_from_treatment_start

        1
        |
        | many
        v
cell_counts
-----------
id             PK
sample_id      FK -> samples.sample_id
population
count
```

`cell_counts` uses a long/tidy representation: one row per sample and cell population. This avoids adding a new column to the schema whenever a new cell population is introduced.

The subject/sample split also avoids repeatedly storing patient-level attributes on every timepoint. Primary keys and indexes on subject and sample foreign keys support joins as the dataset grows.

### Scaling rationale

For hundreds of projects and thousands or millions of samples, the same separation remains useful:

- A project dimension can be split into its own table if project-level metadata grows.
- Subject metadata remains in `subjects`, while repeated measurements remain in `samples` and `cell_counts`.
- Indexes on foreign keys and common analytical filters can support joins and cohort selection.
- Additional analytical results can be stored in separate result tables keyed to a project, analysis run, sample, or subject rather than modifying the core measurement table.
- For substantially larger production workloads, the same logical model can move from SQLite to PostgreSQL without changing the analytical concepts.

## Analysis

### Part 2 — relative frequency

For every sample and each of the five supplied cell populations, the pipeline calculates:

`population count / total count across the five populations × 100`

The result is written to:

`output/part2_frequency_table.csv`

An aggregate composition plot is also written to:

`output/part2_average_composition.html`

### Part 3 — responders vs non-responders

The responder analysis uses melanoma PBMC samples from subjects treated with miraclib and with a recorded response.

Each subject has repeated timepoints. To avoid treating repeated measurements from the same subject as independent observations, the analysis first averages each subject's population frequency across their available samples. It then compares responders and non-responders using a two-sided Mann–Whitney U test, with Benjamini–Hochberg correction across the five populations.

Outputs:

- `output/part3_subject_level_frequencies.csv`
- `output/part3_statistical_results.csv`
- `output/part3_responders_vs_nonresponders.html`

### Part 4 — baseline cohort

The baseline cohort is:

- condition = melanoma
- sample type = PBMC
- treatment = miraclib
- time from treatment start = 0

Outputs:

- `output/part4_baseline_subset.csv`
- `output/part4_by_project.csv`
- `output/part4_by_response.csv`
- `output/part4_by_sex.csv`

## Requested answers

The pipeline writes a machine-readable copy of the answers to `output/requested_answers.csv`.

For clarity, the verified values from the supplied dataset are:

| Question | Answer |
|---|---:|
| Male subjects across the full dataset | **1,810** |
| Female subjects across the full dataset | **1,690** |
| Male subjects in the baseline miraclib/melanoma/PBMC cohort | **344** |
| Female subjects in the baseline miraclib/melanoma/PBMC cohort | **312** |
| Average B-cell count for melanoma males who were responders at time = 0, all sample/treatment types | **10,206.15** |

The last calculation uses the raw `b_cell` count, not the relative-frequency percentage. It filters only melanoma, male, responder, and time 0, leaving sample type and treatment unrestricted.

## Dashboard

Run:

```bash
make dashboard
```

Then open the forwarded Streamlit port 8501 in GitHub Codespaces, or visit:

http://localhost:8501

The dashboard contains:

1. Frequency overview
2. Responder vs non-responder analysis
3. Baseline cohort explorer
4. The requested subject-sex counts and B-cell answer

## Reproducibility

`make pipeline` always rebuilds `cell_count.db` from the checked-in `cell-count.csv`, so a fresh Codespace does not depend on a previously created database.

The generated files under `output/` are reproducible from the input CSV and the Python source.
