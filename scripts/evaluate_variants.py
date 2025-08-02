#!/usr/bin/env python
import vcf

truth_path = snakemake.input["truth_vcf"]
pred_path = snakemake.input["pred_vcf"]
report_path = snakemake.output["benchmark_report"]
qs = snakemake.params["qs"]

def read_variants(path):
    vcf_reader = vcf.Reader(filename=path)
    variant_set = set()
    for rec in vcf_reader:
        # Sadece SNP’leri karşılaştıralım (opsiyonel ama önerilir)
        if len(rec.REF) == 1 and all(len(str(alt)) == 1 for alt in rec.ALT):
            try:
                variant_set.add((rec.CHROM.strip(), rec.POS, rec.REF, str(rec.ALT[0]).strip()))
            except Exception as e:
                print(f"Skipping record at {rec.CHROM}:{rec.POS} due to error: {e}")
    return variant_set

# # Dosya yolları (manuel test için)
# truth_path = "/home/asgurkan/Documents/massive_bioinformatics/variant_calling_benchmark/simulated/simulated_truth.refseq2simseq.SNP.vcf"
# pred_path = "/home/asgurkan/Documents/massive_bioinformatics/variant_calling_benchmark/backup/medaka/medaka_out_good/medaka.vcf.gz"
# report_path = "/home/asgurkan/Documents/massive_bioinformatics/variant_calling_benchmark/benchmark/bad_medaka_metrics_custom.txt"
# qs = "bad"

# Setleri oku
truth_set = read_variants(truth_path)
pred_set = read_variants(pred_path)

# Karşılaştır
tp = truth_set & pred_set
fp = pred_set - truth_set
fn = truth_set - pred_set

# Metrikler
precision = len(tp) / (len(tp) + len(fp)) if (len(tp) + len(fp)) > 0 else 0
recall = len(tp) / (len(tp) + len(fn)) if (len(tp) + len(fn)) > 0 else 0
f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

# Yazdır
report = (
    f"Quality scenario: {qs}\n"
    f"TP: {len(tp)}\n"
    f"FP: {len(fp)}\n"
    f"FN: {len(fn)}\n"
    f"Precision: {precision:.4f}\n"
    f"Recall: {recall:.4f}\n"
    f"F1-score: {f1:.4f}\n"
)

print(report)
with open(report_path, "w") as f:
    f.write(report)
