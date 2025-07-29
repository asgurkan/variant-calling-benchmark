import vcf

QS = "bad"  # kalite etiketi (örnek: bad, medium, good)

def read_variants(path):
    vcf_reader = vcf.Reader(filename=path)
    return {(rec.CHROM, rec.POS, rec.REF, str(rec.ALT[0])) for rec in vcf_reader}

# VCF dosyalarını oku
truth_set = read_variants("simulated_truth.refseq2simseq.SNP.vcf.gz")
longshot_set = read_variants(f"/home/asgurkan/Documents/massive_bioinformatics/variant_calling_benchmark/medaka_out_bad/medaka.vcf.gz")

# Hesapla
tp = truth_set & longshot_set
fp = longshot_set - truth_set
fn = truth_set - longshot_set

# Metrikler
precision = len(tp) / (len(tp) + len(fp)) if (len(tp) + len(fp)) > 0 else 0
recall = len(tp) / (len(tp) + len(fn)) if (len(tp) + len(fn)) > 0 else 0
f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

# Yazdır
report = (
    f"Quality scenario: {QS}\n"
    f"TP: {len(tp)}\n"
    f"FP: {len(fp)}\n"
    f"FN: {len(fn)}\n"
    f"Precision: {precision:.4f}\n"
    f"Recall: {recall:.4f}\n"
    f"F1-score: {f1:.4f}\n"
)

print(report)

# Dosyaya yaz
with open(f"{QS}_metrics_medaka.txt", "w") as f:
    f.write(report)
