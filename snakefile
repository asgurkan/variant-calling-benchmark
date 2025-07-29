

RQ = "bad"
configfile: "config.yaml"

rule all:
    input:
        expand("alignment/{RQ}/aligned_{RQ}.bam", RQ = RQ)

rule simug:
    input: 
        ref_genome = config["ref_genome"]
    output:
        fasta = "simulated/simulated_truth.simseq.genome.fa",
        vcf = "simulated/simulated_truth.refseq2simseq.SNP.vcf",
        mapfile = "simulated/simulated_truth.refseq2simseq.map.txt"
    params:
        prefix = config["simug"]["prefix"],
        simug_script_loc = config["simug"]["simug_script_loc"],
        snp_count = config["simug"]["snp_count"],
        titv_ratio = config["simug"]["titv_ratio"]
    conda:
        "envs/simug.yaml"
    shell:
        """
        perl {params.simug_script_loc} \
            -refseq {input.ref_genome} \
            -snp_count {params.snp_count} \
            -titv_ratio {params.titv_ratio} \
            -prefix {params.prefix}

        """


rule badread:
    input:
        ref_truth = "simulated/simulated_truth.simseq.genome.fa"
    output:
        fastq = "reads/{RQ}_reads.fastq.gz"
    params:
        quantity    = config["badread"][f"{RQ}".strip()]["quantity"],
        error_model = config["badread"][f"{RQ}".strip()]["error_model"],
        qscore_model= config["badread"][f"{RQ}".strip()]["qscore_model"],
        glitches    = config["badread"][f"{RQ}".strip()]["glitches"],
        junk_reads  = config["badread"][f"{RQ}".strip()]["junk_reads"],
        random_reads= config["badread"][f"{RQ}".strip()]["random_reads"],
        chimeras= config["badread"][f"{RQ}".strip()]["chimeras"],
        identity    = config["badread"][f"{RQ}".strip()]["identity"],
        length      = config["badread"][f"{RQ}".strip()]["length"],
        start_adapter_seq = config["badread"][f"{RQ}".strip()]["start_adapter_seq"],
        end_adapter_seq   = config["badread"][f"{RQ}".strip()]["end_adapter_seq"],
        seed        = config["badread"][f"{RQ}".strip()]["seed"]
    conda:
        "envs/badread.yaml"
    shell:
        """
        badread simulate \
            --reference {input.ref_truth} \
            --quantity {params.quantity} \
            --error_model {params.error_model} \
            --qscore_model {params.qscore_model} \
            --glitches {params.glitches} \
            --junk_reads {params.junk_reads} \
            --random_reads {params.random_reads} \
            --chimeras {params.chimeras} \
            --identity {params.identity} \
            --length {params.length} \
            --start_adapter_seq {params.start_adapter_seq} \
            --end_adapter_seq {params.end_adapter_seq} \
            --seed {params.seed} \
            | gzip > reads_good.fastq.gz
        """

rule alignment:
    input:
        reads = "reads/{RQ}_reads.fastq.gz",
        ref = config["ref_genome"]
    output:
        bam = "alignment/{RQ}/aligned_{RQ}.bam"
    params:
        preset = "map-ont"
    conda:
        "envs/alignment.yaml"
    threads: 8
    shell:
        """
        minimap2 -t {threads} -ax {params.preset} {input.ref} {input.reads} \
        | samtools sort -@ {threads} -o {output.bam}

        samtools index {output.bam}
        """
