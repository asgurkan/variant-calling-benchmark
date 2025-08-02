#!/usr/bin/env python
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

benchmark_files = snakemake.input.benchmark_reports
out_dir = snakemake.output.out_dir

os.makedirs(out_dir, exist_ok=True)

data = []
for path in benchmark_files:
    filename = os.path.basename(path)
    tool = "Longshot" if "longshot" in filename.lower() else "Medaka"  # Burada tool belirleniyor
    for_quality = "good" if "good" in filename.lower() else "medium" if "medium" in filename.lower() else "bad"

    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]
        entry = {"File": filename, "Tool": tool, "Quality": for_quality}
        for line in lines:
            if ":" in line:
                key, val = line.split(":")
                key = key.strip()
                try:
                    val = float(val.strip())
                except ValueError:
                    continue
                entry[key] = val
        data.append(entry)

df = pd.DataFrame(data)

# Save CSV summary
df.to_csv(os.path.join(out_dir, "benchmark_metrics_summary.csv"), index=False)

# Clean categories
df["Tool"] = pd.Categorical(df["Tool"], categories=["Longshot", "Medaka"], ordered=True)
df["Quality"] = pd.Categorical(df["Quality"], categories=["good", "medium", "bad"], ordered=True)

# Reshape and plot
metrics_df = df.melt(
    id_vars=["Tool", "Quality"],
    value_vars=["Precision", "Recall", "F1-score"],
    var_name="Metric", value_name="Value"
)

for metric in ["Precision", "Recall", "F1-score"]:
    plt.figure(figsize=(8, 5))
    subset = metrics_df[metrics_df["Metric"] == metric]
    sns.barplot(data=subset, x="Quality", y="Value", hue="Tool", palette="Set2", ci=None)
    plt.title(f"{metric} by Tool and Quality")
    plt.ylim(0, 1.05)
    plt.xlabel("Quality Scenario")
    plt.ylabel(metric)
    plt.legend(title="Tool")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{metric.lower()}_barplot.png"))
    plt.close()
