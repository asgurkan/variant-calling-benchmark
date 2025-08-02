

# RQ = "medium" ["bad", "medium", "good"]
configfile: "config.yaml"

rule all:
    input:
        #expand("variants/{RQ}/medaka_variants.vcf", RQ = RQ)
        # expand("benchmark/{RQ}_{tool}_metrics.txt",
        #        RQ=["bad"],
        #        tool=["longshot", "medaka"])
        expand("benchmark/results", RQ=["medium"], tool=["longshot", "medaka"])

# rule simug:
#     input: 
#         ref_genome = config["ref_genome"]
#     output:
#         fasta = "simulated/simulated_truth.simseq.genome.fa",
#         vcf = "simulated/simulated_truth.refseq2simseq.SNP.vcf.gz",
#         tbi = "simulated/simulated_truth.refseq2simseq.SNP.vcf.gz.tbi",
#         mapfile = "simulated/simulated_truth.refseq2simseq.map.txt"
#     params:
#         prefix = config["simug"]["prefix"],
#         simug_script_loc = config["simug"]["simug_script_loc"],
#         snp_count = config["simug"]["snp_count"],
#         titv_ratio = config["simug"]["titv_ratio"]
#     conda:
#         "envs/simug.yaml"
#     shell:
#         """
#         perl {params.simug_script_loc} \
#             -refseq {input.ref_genome} \
#             -snp_count {params.snp_count} \
#             -titv_ratio {params.titv_ratio} \
#             -prefix {params.prefix}

#         bgzip -f simulated/simulated_truth.refseq2simseq.SNP.vcf
#         tabix -f simulated/simulated_truth.refseq2simseq.SNP.vcf.gz
#         """


# rule badread:
#     input:
#         ref_truth = "simulated/simulated_truth.simseq.genome.fa"
#     output:
#         fastq = "reads/{RQ}_reads.fastq.gz"
#     params:
#         quantity    = config["badread"][f"{RQ}".strip()]["quantity"],
#         error_model = config["badread"][f"{RQ}".strip()]["error_model"],
#         qscore_model= config["badread"][f"{RQ}".strip()]["qscore_model"],
#         glitches    = config["badread"][f"{RQ}".strip()]["glitches"],
#         junk_reads  = config["badread"][f"{RQ}".strip()]["junk_reads"],
#         random_reads= config["badread"][f"{RQ}".strip()]["random_reads"],
#         chimeras    = config["badread"][f"{RQ}".strip()]["chimeras"],
#         identity    = config["badread"][f"{RQ}".strip()]["identity"],
#         length      = config["badread"][f"{RQ}".strip()]["length"],
#         start_adapter_seq = config["badread"][f"{RQ}".strip()]["start_adapter_seq"],
#         end_adapter_seq   = config["badread"][f"{RQ}".strip()]["end_adapter_seq"],
#         seed        = config["badread"][f"{RQ}".strip()]["seed"]
#     conda:
#         "envs/badread.yaml"
#     shell:
#         """
#         badread simulate \
#             --reference {input.ref_truth} \
#             --quantity {params.quantity} \
#             --error_model {params.error_model} \
#             --qscore_model {params.qscore_model} \
#             --glitches {params.glitches} \
#             --junk_reads {params.junk_reads} \
#             --random_reads {params.random_reads} \
#             --chimeras {params.chimeras} \
#             --identity {params.identity} \
#             --length {params.length} \
#             --start_adapter_seq {params.start_adapter_seq} \
#             --end_adapter_seq {params.end_adapter_seq} \
#             --seed {params.seed} \
#             | gzip > {output.fastq}
#         """

# rule alignment:
#     input:
#         reads = "reads/{RQ}_reads.fastq.gz",
#         ref = config["ref_genome"]
#     output:
#         bam = "alignment/{RQ}/aligned_{RQ}.bam"
#     params:
#         preset = "map-ont"
#     conda:
#         "envs/alignment.yaml"
#     threads: 8
#     shell:
#         """
#         minimap2 -t {threads} -ax {params.preset} {input.ref} {input.reads} \
#         | samtools sort -@ {threads} -o {output.bam}

#         samtools index {output.bam}
#         """

# rule medaka_variant_calling:
#     input:
#         bam = "alignment/{RQ}/aligned_{RQ}.bam",
#         ref = config["ref_genome"]
#     output:
#         vcf = "variants/{RQ}/medaka_variants.vcf.gz",
#         tbi = "variants/{RQ}/medaka_variants.vcf.gz.tbi"
#     conda:
#         "envs/medaka.yaml"
#     threads: 12
#     shell:
#         """
#         medaka_variant \
#             -i {input.bam} \
#             -r {input.ref} \
#             -t {threads} \
#             -o variants/{wildcards.RQ}/

#         bcftools sort variants/{wildcards.RQ}/medaka.annotated.vcf -Oz -o {output.vcf}
#         tabix -f {output.vcf}
#         """

# rule longshot_variant_calling:
#     input:
#         bam = "alignment/{RQ}/aligned_{RQ}.bam",
#         ref = config["ref_genome"]
#     output:
#         vcf = "variants/{RQ}/longshot_variants.vcf.gz",
#         tbi = "variants/{RQ}/longshot_variants.vcf.gz.tbi"
#     conda:
#         "envs/longshot.yaml"
#     threads: 8
#     shell:
#         """
#         mkdir -p variants/{wildcards.RQ}

#         longshot \
#             --bam {input.bam} \
#             --ref {input.ref} \
#             --out variants/{wildcards.RQ}/longshot_variants.vcf 


#         bgzip -f variants/{wildcards.RQ}/longshot_variants.vcf
#         tabix -f variants/{wildcards.RQ}/longshot_variants.vcf.gz
#         """

# rule evaluate_variants:
#     input:
#         truth_vcf = "simulated/simulated_truth.refseq2simseq.SNP.vcf.gz",
#         pred_vcf = "variants/{RQ}/{tool}_variants.vcf.gz"
#     output:
#         benchmark_report = "benchmark/{RQ}_{tool}_metrics.txt"
#     params:
#         qs = lambda wildcards: wildcards.RQ
#     conda:
#         "envs/variant_eval.yaml"
#     script:
#         "scripts/evaluate_variants.py"

# rule benchmark_vis:
#     input:
#         benchmark_report = "benchmark/{RQ}_{tool}_metrics.txt"
#     output:
#         out_dir = "benchmark/results/{RQ}_{tool}_plot.png"
#     conda:
#         "envs/variant_eval.yaml"
#     script:
#         "scripts/benchmark_vis.py"

rule benchmark_vis:
    input:
        benchmark_reports = expand("benchmark/{RQ}_{tool}_metrics.txt", RQ=["medium"], tool=["longshot", "medaka"] )
    output:
        out_dir = directory("benchmark/results")
    conda:
        "envs/benchmark_vis.yaml"
    script:
        "scripts/benchmark_vis.py"
