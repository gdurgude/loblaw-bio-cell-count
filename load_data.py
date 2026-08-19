#!/usr/bin/env python3
"""
load_data.py

Initializes a SQLite database (cell_count.db) with a normalized relational
schema and loads all rows from cell-count.csv into it.

Usage:
    python load_data.py

No command-line arguments required. Run this from the repository root,
next to cell-count.csv.
"""

import csv
import os
import sqlite3
import sys

DB_FILE = "cell_count.db"
CSV_FILE = "cell-count.csv"

# The five immune cell population columns present in cell-count.csv.
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

# Alternate header names we tolerate, mapped to the canonical name we use
# internally. This makes the loader resilient to minor naming differences
# (e.g. "gender" vs "sex", "indication" vs "condition").
HEADER_ALIASES = {
    "sex": "sex",
    "gender": "sex",
    "condition": "condition",
    "indication": "condition",
    "disease": "condition",
    "subject": "subject",
    "subject_id": "subject",
    "patient": "subject",
    "patient_id": "subject",
    "project": "project",
    "sample": "sample",
    "sample_id": "sample",
    "sample_type": "sample_type",
    "treatment": "treatment",
    "response": "response",
    "time_from_treatment_start": "time_from_treatment_start",
    "age": "age",
}

SCHEMA = """
CREATE TABLE IF NOT EXISTS subjects (
    subject_id  TEXT PRIMARY KEY,
    project     TEXT,
    condition   TEXT,
    age         INTEGER,
    sex         TEXT,
    treatment   TEXT,
    response    TEXT
);

CREATE TABLE IF NOT EXISTS samples (
    sample_id                  TEXT PRIMARY KEY,
    subject_id                 TEXT NOT NULL,
    sample_type                TEXT,
    time_from_treatment_start  REAL,
    FOREIGN KEY (subject_id) REFERENCES subjects (subject_id)
);

CREATE TABLE IF NOT EXISTS cell_counts (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id   TEXT NOT NULL,
    population  TEXT NOT NULL,
    count       INTEGER NOT NULL,
    FOREIGN KEY (sample_id) REFERENCES samples (sample_id)
);

CREATE INDEX IF NOT EXISTS idx_samples_subject ON samples (subject_id);
CREATE INDEX IF NOT EXISTS idx_cellcounts_sample ON cell_counts (sample_id);
"""


def normalize_headers(fieldnames):
    """Map raw CSV headers (case/spacing-insensitive) to canonical names."""
    mapping = {}
    for raw in fieldnames:
        key = raw.strip().lower().replace(" ", "_")
        if key in POPULATIONS:
            mapping[raw] = key
        elif key in HEADER_ALIASES:
            mapping[raw] = HEADER_ALIASES[key]
        else:
            mapping[raw] = key
    return mapping


def init_db(conn):
    conn.executescript(SCHEMA)
    conn.commit()


def load_csv(conn, csv_path):
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"{csv_path} appears to be empty or has no header row.")

        header_map = normalize_headers(reader.fieldnames)

        if "sample" not in header_map.values():
            raise ValueError(
                f"{csv_path} has no 'sample' column (found: {reader.fieldnames}). "
                "Cannot load without a sample identifier."
            )
        missing_populations = [p for p in POPULATIONS if p not in header_map.values()]
        if missing_populations:
            raise ValueError(
                f"{csv_path} is missing expected population column(s): {missing_populations}"
            )

        subjects_seen = set()
        n_samples = 0
        n_cell_rows = 0
        n_skipped = 0

        cur = conn.cursor()

        for line_num, raw_row in enumerate(reader, start=2):  # header is line 1
            row = {header_map[k]: v for k, v in raw_row.items() if k in header_map}

            sample_id = (row.get("sample") or "").strip()
            if not sample_id:
                print(f"WARNING: skipping line {line_num}, no sample id", file=sys.stderr)
                n_skipped += 1
                continue

            try:
                # A subject id is expected in most versions of this dataset. If
                # the file doesn't provide one, fall back to using the sample id
                # itself (i.e. treat each sample as its own subject).
                subject_id = (row.get("subject") or "").strip() or sample_id

                if subject_id not in subjects_seen:
                    age_val = row.get("age")
                    age_val = int(age_val) if age_val not in (None, "") else None
                    cur.execute(
                        """
                        INSERT OR IGNORE INTO subjects
                            (subject_id, project, condition, age, sex, treatment, response)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            subject_id,
                            row.get("project"),
                            row.get("condition"),
                            age_val,
                            row.get("sex"),
                            row.get("treatment"),
                            (row.get("response") or "").strip().lower() or None,
                        ),
                    )
                    subjects_seen.add(subject_id)

                tfts = row.get("time_from_treatment_start")
                tfts_val = float(tfts) if tfts not in (None, "") else None

                cur.execute(
                    """
                    INSERT OR REPLACE INTO samples
                        (sample_id, subject_id, sample_type, time_from_treatment_start)
                    VALUES (?, ?, ?, ?)
                    """,
                    (sample_id, subject_id, row.get("sample_type"), tfts_val),
                )
                n_samples += 1

                for pop in POPULATIONS:
                    if pop in row and row[pop] not in (None, ""):
                        cur.execute(
                            """
                            INSERT INTO cell_counts (sample_id, population, count)
                            VALUES (?, ?, ?)
                            """,
                            (sample_id, pop, int(float(row[pop]))),
                        )
                        n_cell_rows += 1

            except (ValueError, TypeError) as e:
                print(f"WARNING: skipping line {line_num} ({sample_id}): {e}", file=sys.stderr)
                n_skipped += 1
                continue

        conn.commit()
        if n_skipped:
            print(f"Skipped {n_skipped} malformed row(s), see warnings above.", file=sys.stderr)
        return len(subjects_seen), n_samples, n_cell_rows


def main():
    if not os.path.exists(CSV_FILE):
        print(f"ERROR: {CSV_FILE} not found in the current directory.", file=sys.stderr)
        sys.exit(1)

    # Start fresh each run so re-loading is idempotent.
    if os.path.exists(DB_FILE):
        os.remove(DB_FILE)

    conn = sqlite3.connect(DB_FILE)
    try:
        init_db(conn)
        n_subjects, n_samples, n_cell_rows = load_csv(conn, CSV_FILE)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        conn.close()
        if os.path.exists(DB_FILE):
            os.remove(DB_FILE)
        sys.exit(1)
    else:
        print(f"Loaded {n_subjects} subjects, {n_samples} samples, "
              f"{n_cell_rows} cell count records into {DB_FILE}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
