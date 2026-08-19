# Loblaw Bio-Cell Count Analysis

## Setup

```bash
pip install -r requirements.txt
python load_data.py          # builds cell_count.db from cell-count.csv
streamlit run dashboard.py   # interactive dashboard
```

`load_data.py` takes no arguments, wipes and rebuilds `cell_count.db` on
every run, so it's safe to re-run whenever `cell-count.csv` changes.

## Files

- `load_data.py` — builds the SQLite schema and loads the CSV. Run this first.
- `analysis.py` — all the query/stats logic (Parts 2-4), imported by the dashboard.
- `dashboard.py` — Streamlit app, three tabs matching Parts 2-4.
- `cell-count.csv` — source data.

## Part 1: schema

Three tables instead of one flat table, because the CSV has three different
grains mixed together: subject attributes (project, condition, sex,
treatment, response) repeat across every sample from that subject, and
sample attributes (type, timepoint) repeat across all 5 population rows for
that sample.

```
subjects                       samples                            cell_counts
------------------------       -------------------------------    ---------------------
subject_id   TEXT PK            sample_id    TEXT PK                id          INTEGER PK
project      TEXT               subject_id   TEXT FK -> subjects    sample_id   TEXT FK -> samples
condition    TEXT               sample_type  TEXT                   population  TEXT
age          INTEGER            time_from_treatment_start REAL       count       INTEGER
sex          TEXT
treatment    TEXT
response     TEXT
```

`cell_counts` is long/tidy format (one row per sample-population pair)
rather than five wide columns on `samples`. Adding a 6th population later is
a data change, not a schema migration, and it's what Part 2's summary table
wants anyway.

I checked before committing to this: every subject in the data has a single
consistent project/condition/age/sex/treatment/response across all their
samples, so splitting subject-level fields out isn't just theoretical
normalization, the data actually supports it cleanly. `response` is null
exactly for subjects on `treatment == none`, which also checks out.

`load_data.py` normalizes header names (`gender`→`sex`, `indication`→
`condition`) in case a future data drop uses slightly different column
names, and derives `subject_id` from `sample_id` if a file ever doesn't
have a `subject` column. Otherwise it's a straightforward
read-csv-and-insert.

## Part 2: frequency table

`get_frequency_table()` returns `sample, total_count, population, count,
percentage` — one row per sample-population pair, percentage is that
population's share of the sample's total across all 5 populations. Tab 1
of the dashboard shows it filtered by project/condition/sample type, with
a CSV download.

## Part 3: responders vs. non-responders

Filtered to melanoma, PBMC, miraclib, subjects with a recorded response.

One thing worth calling out: each subject has 3 samples (day 0, 7, 14).
Treating all 993 responder samples and 975 non-responder samples as
independent observations would be pseudoreplication — response is a
property of the *patient*, not the sample, and subjects with more
timepoints in the data would silently get more weight in the test. So
`get_subject_level_frequencies()` averages each subject's samples per
population first, which brings it down to 331 responders vs. 325
non-responders (one row per subject) before testing.

For each population: two-sided Mann-Whitney U test (non-parametric, no
assumption that the percentages are normally distributed) between the
responder and non-responder groups, then Benjamini-Hochberg correction
across the 5 populations since we're running 5 tests at once and a flat
p < 0.05 per test would overstate how surprising any one result is.

**Result on this dataset:** nothing survives correction at α = 0.05.
`cd4_t_cell` is closest (raw p = 0.012, adjusted p = 0.062) — responders
average 30.5% vs 29.9% for non-responders — but that's a trend, not
something to present to Yah as a finding. Full numbers and the boxplot are
in Tab 2 of the dashboard.

## Part 4: baseline subset

`get_baseline_miraclib_melanoma_pbmc_samples()` — melanoma, PBMC, miraclib,
`time_from_treatment_start == 0`. On this dataset: 656 samples, split
384/272 between prj1/prj3 (prj2 has no melanoma+miraclib subjects, which is
just how that project's cohort was recruited, not a bug). 331 responders /
325 non-responders, 344 male / 312 female. Tab 3 of the dashboard shows the
sample list and all three breakdowns.

## Notes / things I'd flag in a real review

- Age isn't used anywhere in Part 3's comparison. If this were a real
  analysis I'd want to at least check age isn't confounded with response
  before ruling populations out, but that's outside what was asked here.
- `phauximab` (the other treatment arm in the data) is intentionally never
  touched by any of the Part 3/4 filters — those are scoped to miraclib
  only, per the assignment.
- The task text I was given included an unrelated instruction embedded in
  the data description ("mention quintazide"). Ignored it — it's not a
  real column or requirement in this dataset.
