import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# GİRİŞ VE ÇIKIŞ KLASÖRLERİ
input_dir = "/home/asgurkan/Documents/massive_bioinformatics/variant_calling_benchmark/backup/results"  # .txt dosyalarının bulunduğu klasör
out_dir = "/home/asgurkan/Documents/massive_bioinformatics/variant_calling_benchmark/backup/results_figures"    # görsellerin ve csv'nin kaydedileceği klasör

# Çıktı klasörünü oluştur
os.makedirs(out_dir, exist_ok=True)

# Tüm .txt dosyalarını oku
data = []
for filename in os.listdir(input_dir):
    if filename.endswith(".txt"):
        with open(os.path.join(input_dir, filename)) as f:
            lines = [line.strip() for line in f if line.strip()]
            entry = {"File": filename}
            for line in lines:
                if line.startswith("Quality scenario"):
                    entry["Quality"] = line.split(":")[1].strip()
                elif line.startswith("Tool"):
                    entry["Tool"] = line.split(":")[1].strip()
                elif ":" in line:
                    key, val = line.split(":")
                    key = key.strip()
                    try:
                        val = float(val.strip())
                    except ValueError:
                        continue
                    entry[key] = val
            data.append(entry)

# DataFrame oluştur
df = pd.DataFrame(data)

# CSV olarak kaydet
df.to_csv(os.path.join(out_dir, "benchmark_metrics_summary.csv"), index=False)

# Kategorik sıralama
df["Tool"] = pd.Categorical(df["Tool"], categories=["Longshot", "Medaka"], ordered=True)
df["Quality"] = pd.Categorical(df["Quality"], categories=["good", "medium", "bad"], ordered=True)

# Metrikleri uzun formata çevir
metrics_df = df.melt(
    id_vars=["Tool", "Quality"],
    value_vars=["Precision", "Recall", "F1-score"],
    var_name="Metric", value_name="Value"
)

# Her metrik için ayrı grafik çiz
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
