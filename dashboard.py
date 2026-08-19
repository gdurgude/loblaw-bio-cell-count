"""
dashboard.py

Interactive Streamlit dashboard for Bob Loblaw's immune cell population
analysis. Run with:

    streamlit run dashboard.py

Requires cell_count.db to already exist (run `python load_data.py` first).
"""

import os

import pandas as pd
import plotly.express as px
import streamlit as st

import analysis as A

st.set_page_config(page_title="Loblaw Bio | Cell Count Analysis", layout="wide")

DB_FILE = "cell_count.db"

st.title("🧬 Loblaw Bio — Immune Cell Population Dashboard")
st.caption(
    "Cell count analysis for miraclib clinical trial samples. "
    "Data source: cell-count.csv, loaded via load_data.py."
)

if not os.path.exists(DB_FILE):
    st.error(
        f"Database file `{DB_FILE}` not found. Run `python load_data.py` "
        "in the repository root first, then reload this page."
    )
    st.stop()

conn = A.get_connection(DB_FILE)

tab1, tab2, tab3 = st.tabs([
    "1. Frequency Overview",
    "2. Responders vs Non-Responders",
    "3. Baseline Subset Explorer",
])

# ---------------------------------------------------------------------------
# Tab 1: Part 2 — frequency table
# ---------------------------------------------------------------------------
with tab1:
    st.header("Relative frequency of each cell population, per sample")
    freq = A.get_frequency_table(conn)
    meta = A.get_full_sample_table(conn)

    with st.expander("Filter samples", expanded=False):
        col1, col2, col3 = st.columns(3)
        projects = sorted(meta["project"].dropna().unique())
        conditions = sorted(meta["condition"].dropna().unique())
        sample_types = sorted(meta["sample_type"].dropna().unique())

        sel_projects = col1.multiselect("Project", projects, default=projects)
        sel_conditions = col2.multiselect("Condition", conditions, default=conditions)
        sel_types = col3.multiselect("Sample type", sample_types, default=sample_types)

    filtered_samples = meta[
        meta["project"].isin(sel_projects)
        & meta["condition"].isin(sel_conditions)
        & meta["sample_type"].isin(sel_types)
    ]["sample_id"]

    display_freq = freq[freq["sample"].isin(filtered_samples)]

    st.dataframe(display_freq, use_container_width=True, hide_index=True)
    st.download_button(
        "Download table as CSV",
        display_freq.to_csv(index=False),
        file_name="cell_frequencies.csv",
        mime="text/csv",
    )

    st.subheader("Average population composition (filtered samples)")
    avg_comp = display_freq.groupby("population")["percentage"].mean().reset_index()
    fig = px.bar(
        avg_comp, x="population", y="percentage",
        labels={"percentage": "Mean % of total cells", "population": "Population"},
    )
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("Look up a single sample")
    sample_choice = st.selectbox("Sample", sorted(display_freq["sample"].unique()))
    sample_detail = display_freq[display_freq["sample"] == sample_choice]
    fig2 = px.pie(sample_detail, names="population", values="percentage",
                  title=f"{sample_choice} composition")
    st.plotly_chart(fig2, use_container_width=True)

# ---------------------------------------------------------------------------
# Tab 2: Part 3 — statistical analysis
# ---------------------------------------------------------------------------
with tab2:
    st.header("Melanoma, PBMC, miraclib: responders vs non-responders")
    st.markdown(
        "Comparing relative frequency of each immune cell population between "
        "**responders** and **non-responders**, restricted to melanoma PBMC "
        "samples from subjects treated with miraclib."
    )

    data = A.get_subject_level_frequencies(conn)

    if data.empty:
        st.warning(
            "No melanoma / PBMC / miraclib samples with a recorded response "
            "were found in this dataset."
        )
    else:
        st.subheader("Boxplots by population")
        st.caption(
            "One point per subject (samples averaged across timepoints), so "
            "subjects aren't double-counted."
        )
        fig = px.box(
            data,
            x="population",
            y="percentage",
            color="response",
            points="all",
            category_orders={"response": ["yes", "no"]},
            labels={"percentage": "Relative frequency (%)", "population": "Population",
                    "response": "Response"},
            color_discrete_map={"yes": "#2E86AB", "no": "#C73E1D"},
        )
        st.plotly_chart(fig, use_container_width=True)

        st.subheader("Statistical test results (Mann-Whitney U, BH-adjusted)")
        stat_table = A.compare_responders_vs_nonresponders(conn)
        st.dataframe(stat_table, use_container_width=True, hide_index=True)

        sig = stat_table[stat_table["significant"] == True]
        if len(sig):
            pops = ", ".join(sig["population"].tolist())
            st.success(f"**Significant differences (FDR-adjusted p < 0.05):** {pops}")
        else:
            top = stat_table.iloc[0]
            st.info(
                f"No population reaches significance after FDR correction. Closest is "
                f"**{top['population']}** (raw p={top['p_value']:.4f}, adjusted "
                f"p={top['p_adj']:.4f}) — a trend worth watching as more data comes in, "
                "but not something to bring to Yah as a finding yet."
            )

        st.caption(
            "Test: two-sided Mann-Whitney U test per population (non-parametric, "
            "appropriate given small clinical sample sizes and no normality "
            "assumption). P-values adjusted across the 5 populations using the "
            "Benjamini-Hochberg FDR procedure."
        )

# ---------------------------------------------------------------------------
# Tab 3: Part 4 — baseline subset explorer
# ---------------------------------------------------------------------------
with tab3:
    st.header("Baseline miraclib melanoma PBMC samples")
    st.markdown(
        "Melanoma PBMC samples at baseline (`time_from_treatment_start == 0`) "
        "from subjects treated with miraclib."
    )

    summ = A.summarize_baseline_subset(conn)
    subset = summ["subset"]

    st.metric("Matching samples", len(subset))

    st.dataframe(
        subset[["sample_id", "subject_id", "project", "sex", "response",
                "time_from_treatment_start"]],
        use_container_width=True,
        hide_index=True,
    )

    col1, col2, col3 = st.columns(3)
    with col1:
        st.subheader("By project")
        st.dataframe(summ["by_project"], hide_index=True, use_container_width=True)
        if len(summ["by_project"]):
            st.plotly_chart(
                px.bar(summ["by_project"], x="project", y="n_samples"),
                use_container_width=True,
            )
    with col2:
        st.subheader("By response")
        st.dataframe(summ["by_response"], hide_index=True, use_container_width=True)
        if len(summ["by_response"]):
            st.plotly_chart(
                px.pie(summ["by_response"], names="response", values="n_subjects"),
                use_container_width=True,
            )
    with col3:
        st.subheader("By sex")
        st.dataframe(summ["by_sex"], hide_index=True, use_container_width=True)
        if len(summ["by_sex"]):
            st.plotly_chart(
                px.pie(summ["by_sex"], names="sex", values="n_subjects"),
                use_container_width=True,
            )

conn.close()
