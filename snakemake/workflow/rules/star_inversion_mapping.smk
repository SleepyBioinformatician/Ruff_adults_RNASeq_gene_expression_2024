# Author: Alex Zemella
# Copyright: Copyright 2023, Alex Zemella
# Email: alex.zemella@bi.mpg.de
# License: Max Planck Institute for Biological Intelligence (MPIBI)

##################################################################################################################################################
### STAR mapping of Faeder and Satellite RNA-Seq reads on each inversion morph-specific consensus genome
### In step 1 I compressed Faeder consensus INDELs VCF
### In step 2 I indexed Faeder consensus INDELs VCF
### In step 3 I compressed Faeder consensus SNPs VCF
### In step 4 I indexed Faeder consensus SNPs VCF
### In step 5 I combined Faeder INDELs and SNPs into a single VCF
### In step 6 I compressed Satellite consensus INDELs VCF
### In step 7 I indexed Satellite consensus INDELs VCF
### In step 8 I compressed Satellite consensus SNPs VCF
### In step 9 I indexed Satellite consensus SNPs VCF
### In step 10 I combined Satellite INDELs and SNPs into a single VCF
### In step 11 I mapped Faeder RNA-Seq libraries on the new Faeder's consensus genome using STAR 2nd pass mapping strategy
### In step 12 I mapped Satellite RNA-Seq libraries on the new Satellite's consensus genome using STAR 2nd pass mapping strategy
### The final outputs of this pipeline are alignment BAM files and gene count files for Faeders and Satellites
##################################################################################################################################################

### Step 1
#rule bcftools_bgzip_inversion_hetINDELs_Faeder:
#    input:
#        indels = "gatk4/ruff_adults_variants/consensus_inversion_hetINDELs_faeder.vcf",
#    output:
#        bgzip = "star_consensus/variants/consensus_inversion_hetINDELs_faeder.vcf.gz",
#    threads: 1
#    conda:
#        "../envs/bcftools.yaml",
#    shell:
#        "bgzip -c {input.indels} > {output.bgzip}"

### Step 2
#rule bcftools_tabix_inversion_hetINDELs_Faeder:
#    input:
#        bgzip = "star_consensus/variants/consensus_inversion_hetINDELs_faeder.vcf.gz",
#    output:
#        tabix = "star_consensus/variants/consensus_inversion_hetINDELs_faeder.vcf.gz.tbi",
#    threads: 1
#    priority: 1000
#    conda:
#        "../envs/bcftools.yaml",
#    shell:
#        "tabix -p vcf {input.bgzip}"

### Step 3
#rule bcftools_bgzip_inversion_hetSNPs_Faeder:
#    input:
#        snps = "gatk4/ruff_adults_variants/consensus_inversion_hetSNPs_faeder.vcf",
#    output:
#        bgzip = "star_consensus/variants/consensus_inversion_hetSNPs_faeder.vcf.gz",
#    threads: 1
#    conda:
#        "../envs/bcftools.yaml",
#    shell:
#        "bgzip -c {input.snps} > {output.bgzip}"

### Step 4
#rule bcftools_tabix_inversion_hetSNPs_Faeder:
#    input:
#        bgzip = "star_consensus/variants/consensus_inversion_hetSNPs_faeder.vcf.gz",
#    output:
#        tabix = "star_consensus/variants/consensus_inversion_hetSNPs_faeder.vcf.gz.tbi",
#    threads: 1
#    priority: 1000
#    conda:
#        "../envs/bcftools.yaml",
#    shell:
#        "tabix -p vcf {input.bgzip}"

### Step 5
#rule bcftools_concat_hetVARIANTS_Faeder:
#    input:
#        snps = "star_consensus/variants/consensus_inversion_hetSNPs_faeder.vcf.gz",
#        indels = "star_consensus/variants/consensus_inversion_hetINDELs_faeder.vcf.gz",
#    output:
#        merged = "star_consensus/variants/consensus_inversion_hetSNPs_and_hetINDELs_merged_faeder.vcf",
#    params:
#        extra = "--allow-overlaps",
#    threads: 1
#    conda:
#        "../envs/bcftools.yaml",
#    shell:
#        "bcftools concat {params.extra} {input.snps} {input.indels} -o {output.merged}"

### Step 6
#rule bcftools_bgzip_inversion_hetINDELs_Satellite:
#    input:
#        indels = "gatk4/ruff_adults_variants/consensus_inversion_hetINDELs_satellite.vcf",
#    output:
#        bgzip = "star_consensus/variants/consensus_inversion_hetINDELs_satellite.vcf.gz",
#    threads: 1
#    conda:
#        "../envs/bcftools.yaml",
#    shell:
#        "bgzip -c {input.indels} > {output.bgzip}"

### Step 7
#rule bcftools_tabix_inversion_hetINDELs_Satellite:
#    input:
#        bgzip = "star_consensus/variants/consensus_inversion_hetINDELs_satellite.vcf.gz",
#    output:
#        tabix = "star_consensus/variants/consensus_inversion_hetINDELs_satellite.vcf.gz.tbi",
#    threads: 1
#    priority: 1000
#    conda:
#        "../envs/bcftools.yaml",
#    shell:
#        "tabix -p vcf {input.bgzip}"

### Step 8
#rule bcftools_bgzip_inversion_hetSNPs_Satellite:
#    input:
#        snps = "gatk4/ruff_adults_variants/consensus_inversion_hetSNPs_satellite.vcf",
#    output:
#        bgzip = "star_consensus/variants/consensus_inversion_hetSNPs_satellite.vcf.gz",
#    threads: 1
#    conda:
#        "../envs/bcftools.yaml",
#    shell:
#        "bgzip -c {input.snps} > {output.bgzip}"

### Step 9
#rule bcftools_tabix_inversion_hetSNPs_Satellite:
#    input:
#        bgzip = "star_consensus/variants/consensus_inversion_hetSNPs_satellite.vcf.gz",
#    output:
#        tabix = "star_consensus/variants/consensus_inversion_hetSNPs_satellite.vcf.gz.tbi",
#    threads: 1
#    priority: 1000
#    conda:
#        "../envs/bcftools.yaml",
#    shell:
#        "tabix -p vcf {input.bgzip}"

### Step 10
#rule bcftools_concat_hetVARIANTS_Satellite:
#    input:
#        snps = "star_consensus/variants/consensus_inversion_hetSNPs_satellite.vcf.gz",
#        indels = "star_consensus/variants/consensus_inversion_hetINDELs_satellite.vcf.gz",
#    output:
#        merged = "star_consensus/variants/consensus_inversion_hetSNPs_and_hetINDELs_merged_satellite.vcf",
#    params:
#        extra = "--allow-overlaps",
#    threads: 1
#    conda:
#        "../envs/bcftools.yaml",
#    shell:
#        "bcftools concat {params.extra} {input.snps} {input.indels} -o {output.merged}"

### Step 11
rule STAR_genome_index_Faeder:
    input:
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        annotation = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.gtf",
    output:
        directory("star_consensus/faeder/GenomeDir"),
    params:
        extra1 = "--runMode genomeGenerate",
        extra2 = "--genomeDir star_consensus/faeder/GenomeDir",
        extra3 = "--genomeFastaFiles genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        extra4 = "--sjdbGTFfile genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.gtf",
        extra5 = "--sjdbOverhang 99",
        extra6 = "--runThreadN 10",
        extra7 = "--genomeTransformVCF gatk4/ruff_adults_variants/consensus_inversion_hetSNPs_faeder.vcf",
        extra8 = "--genomeTransformType Haploid",
        extra9 = "--genomeTransformOutput SAM SJ",
    threads: 20
    conda:
        "../envs/star.yaml",
    shell:
        "STAR {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5}"
        " {params.extra6} {params.extra7} {params.extra8} {params.extra9}"

rule STAR_mapping_1Pass_Faeder:
    input:
        fastq1 = "trimgalore/{fae_sample}_R1_val_1.fq.gz",
        fastq2 = "trimgalore/{fae_sample}_R2_val_2.fq.gz",
        index = "star_consensus/faeder/GenomeDir",
    output:
        junctions = "star_consensus/faeder/junctions/{fae_sample}_SJ.out.tab",
    params:
        extra1 = "--runMode alignReads",
        extra2 = "--readFilesCommand zcat",
        extra3 = "--outSAMtype None",
        extra4 = "--runThreadN 5",
        extra5 = "--outFileNamePrefix star_consensus/faeder/junctions/{fae_sample}_",
    threads: 20
    conda:
        "../envs/star.yaml",
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} "
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5}"

rule junctions_filtering_Faeder:
    input:
        SJ = expand("star_consensus/faeder/junctions/{fae_sample}_SJ.out.tab", fae_sample = config["faeder_samples"]),
    output:
        filteredSJ = "star_consensus/faeder/junctions/SJ_faeder_merged_filtered_junctions_inversion.tab",
    threads: 1
    shell:
        "cat {input.SJ} | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > {output.filteredSJ}"

rule STAR_mapping_2Pass_Faeder:
    input:
        fastq1 = "trimgalore/{fae_sample}_R1_val_1.fq.gz",
        fastq2 = "trimgalore/{fae_sample}_R2_val_2.fq.gz",
        index = "star_consensus/faeder/GenomeDir",
        junctions = "star_consensus/faeder/junctions/SJ_faeder_merged_filtered_junctions_inversion.tab",
    output:
        bam = "star_consensus/faeder/mapped/{fae_sample}_Aligned.out.bam",
        counts = "star_consensus/faeder/mapped/{fae_sample}_ReadsPerGene.out.tab",
    params:
        extra1 = "--runMode alignReads",
        extra2 = "--quantMode GeneCounts",
        extra3 = "--readFilesCommand zcat",
        extra4 = "--outSAMtype BAM Unsorted",
        extra5 = "--outSAMmapqUnique 60",
        extra6 = "--outSAMattributes NH HI AS nM NM MD jM jI XS MC ch",
        extra7 = "--runThreadN 5",
        extra8 = "--outSAMattrRGline ID:{fae_sample}",
        extra9 = "--outFileNamePrefix star_consensus/faeder/mapped/{fae_sample}_",
    threads: 20
    conda:
        "../envs/star.yaml",
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} --sjdbFileChrStartEnd {input.junctions}"
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5} {params.extra6}"
        " {params.extra7} {params.extra8} {params.extra9}"

### Step 12
rule STAR_genome_index_Satellite:
    input:
        genome = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        annotation = "genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.gtf",
    output:
        directory("star_consensus/satellite/GenomeDir"),
    params:
        extra1 = "--runMode genomeGenerate",
        extra2 = "--genomeDir star_consensus/satellite/GenomeDir",
        extra3 = "--genomeFastaFiles genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.fna",
        extra4 = "--sjdbGTFfile genome/ncbi_genome/GCF_001431845.1_ASM143184v1_genomic.gtf",
        extra5 = "--sjdbOverhang 99",
        extra6 = "--runThreadN 10",
        extra7 = "--genomeTransformVCF gatk4/ruff_adults_variants/consensus_inversion_hetSNPs_faeder.vcf",
        extra8 = "--genomeTransformType Haploid",
        extra9 = "--genomeTransformOutput SAM SJ",
    threads: 20
    conda:
        "../envs/star.yaml",
    shell:
        "STAR {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5}"
        " {params.extra6} {params.extra7} {params.extra8} {params.extra9}"

rule STAR_mapping_1Pass_Satellite:
    input:
        fastq1 = "trimgalore/{sat_sample}_R1_val_1.fq.gz",
        fastq2 = "trimgalore/{sat_sample}_R2_val_2.fq.gz",
        index = "star_consensus/satellite/GenomeDir",
    output:
        junctions = "star_consensus/satellite/junctions/{sat_sample}_SJ.out.tab",
    params:
        extra1 = "--runMode alignReads",
        extra2 = "--readFilesCommand zcat",
        extra3 = "--outSAMtype None",
        extra4 = "--runThreadN 5",
        extra5 = "--outFileNamePrefix star_consensus/satellite/junctions/{sat_sample}_",
    threads: 20
    conda:
        "../envs/star.yaml",
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} "
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5}"

rule junctions_filtering_Satellite:
    input:
        SJ = expand("star_consensus/satellite/junctions/{sat_sample}_SJ.out.tab", sat_sample = config["satellite_samples"]),
    output:
        filteredSJ = "star_consensus/satellite/junctions/SJ_satellite_merged_filtered_junctions_inversion.tab",
    threads: 1
    shell:
        "cat {input.SJ} | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > {output.filteredSJ}"

rule STAR_mapping_2Pass_Satellite:
    input:
        fastq1 = "trimgalore/{sat_sample}_R1_val_1.fq.gz",
        fastq2 = "trimgalore/{sat_sample}_R2_val_2.fq.gz",
        index = "star_consensus/satellite/GenomeDir",
        junctions = "star_consensus/satellite/junctions/SJ_satellite_merged_filtered_junctions_inversion.tab",
    output:
        bam = "star_consensus/satellite/mapped/{sat_sample}_Aligned.out.bam",
        counts = "star_consensus/satellite/mapped/{sat_sample}_ReadsPerGene.out.tab",
    params:
        extra1 = "--runMode alignReads",
        extra2 = "--quantMode GeneCounts",
        extra3 = "--readFilesCommand zcat",
        extra4 = "--outSAMtype BAM Unsorted",
        extra5 = "--outSAMmapqUnique 60",
        extra6 = "--outSAMattributes NH HI AS nM NM MD jM jI XS MC ch",
        extra7 = "--runThreadN 5",
        extra8 = "--outSAMattrRGline ID:{sat_sample}",
        extra9 = "--outFileNamePrefix star_consensus/satellite/mapped/{sat_sample}_",
    threads: 20
    conda:
        "../envs/star.yaml",
    shell:
        "STAR --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} --sjdbFileChrStartEnd {input.junctions}"
        " {params.extra1} {params.extra2} {params.extra3} {params.extra4} {params.extra5} {params.extra6}"
        " {params.extra7} {params.extra8} {params.extra9}"


