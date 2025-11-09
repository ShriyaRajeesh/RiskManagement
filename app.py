import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import os

st.set_page_config(page_title="Parallel Financial Stress Test Dashboard", layout="wide")

st.title("💹 Parallel Stress Testing Dashboard")
st.markdown("### Simulated Portfolio Risk under Market Stress Scenarios")

# Load summary
if not os.path.exists("scenario_summary.csv"):
    st.warning("Please run the C++ simulation first (stress_test.exe) to generate data.")
else:
    summary = pd.read_csv("scenario_summary.csv")
    st.dataframe(summary.style.highlight_max(axis=0))

    # ---- Bar chart of VaR and CVaR ----
    st.markdown("## 📊 Value at Risk (VaR) and Conditional VaR (CVaR) Comparison")
    import numpy as np

    x = np.arange(len(summary["scenario"]))  # positions
    bar_width = 0.35

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(x - bar_width/2, summary["VaR95"], width=bar_width, label="VaR 95%", alpha=0.7)
    ax.bar(x + bar_width/2, summary["CVaR95"], width=bar_width, label="CVaR 95%", alpha=0.7)

    ax.set_xticks(x)
    ax.set_xticklabels(summary["scenario"])
    ax.set_ylabel("Risk Metric (%)")
    ax.set_title("VaR & CVaR Comparison Across Scenarios")
    ax.legend()
    st.pyplot(fig)

    # ---- Select scenario for histogram ----
    st.markdown("## 🔍 Scenario Return Distribution")
    scenarios = summary["scenario"].tolist()
    selected_scenario = st.selectbox("Select a scenario:", scenarios)

    file_name = f"returns_{selected_scenario}.csv"
    if os.path.exists(file_name):
        returns_df = pd.read_csv(file_name)
        fig2, ax2 = plt.subplots(figsize=(8, 4))
        ax2.hist(returns_df["return"], bins=50, color="skyblue", edgecolor="black")
        ax2.axvline(x=returns_df["return"].quantile(0.05), color="red", linestyle="--", label="VaR 95%")
        ax2.set_title(f"Return Distribution: {selected_scenario}")
        ax2.set_xlabel("Portfolio Return (%)")
        ax2.set_ylabel("Frequency")
        ax2.legend()
        st.pyplot(fig2)
    else:
        st.warning(f"No data found for {selected_scenario}. Run simulation first.")
