import pandas as pd

configfile: "config.yaml" # loads contents of config.yaml file into a global dict called "config"
# required fields in config.yaml:
# - samples: tab-delimited file with columns 'sample','run','fastq'
# - ref_genome: path to (indexed) reference genome
# - adapters: path to .fa with adapter sequences
# optional fields:
# - fastq_dir: common path to fastq files (gets added to all paths defined in the samples fastq column)

DATASETS = pd.read_table(config["samples"], dtype=str).set_index(["sample", "run"])

if config["fastq_dir"]:
    DATASETS['fastq'] = config["fastq_dir"] + DATASETS['fastq'].astype('str')

SAMPLE_LIST = set([c[0] for c in DATASETS.index])
SAMPLE_RUN_COMBS = DATASETS.index.tolist() # list of (sample, run) tuples

print(DATASETS)


def get_fastq(wildcards):
    ''' returns the 2 fastqs for a given sample-run combination'''
    fastq = DATASETS.loc[(wildcards.sample, wildcards.run), ["fastq"]]
    return {"r1": fastq.fastq.replace('*','1'), "r2": fastq.fastq.replace('*','2')}


def get_sample_runs(wildcards):
    '''generates all sample-run combinations for a given sample'''
    runs = list(DATASETS.loc[wildcards.sample].index)
    return expand("filtered/{sample}-{run}.filtered.bam", run=runs, **wildcards)

rule all:
    input:
        ["dedup/{}.dedup.bam".format(s) for s in SAMPLE_LIST]
        #["filtered/{0}-{1}.filtered.bam".format(c[0],c[1]) for c in SAMPLE_RUN_COMBS],
        #["qc/fastqc_raw/{0}-{1}_R1_fastqc.html".format(c[0],c[1]) for c in SAMPLE_RUN_COMBS],
        #["qc/fastqc_trimmed/{0}-{1}.1P_fastqc.html".format(c[0],c[1]) for c in SAMPLE_RUN_COMBS],
        #["qc/mapped_stats/{0}-{1}.txt".format(c[0],c[1]) for c in SAMPLE_RUN_COMBS],
        #["qc/filtered_stats/{0}-{1}.filtered.txt".format(c[0],c[1]) for c in SAMPLE_RUN_COMBS],
        #"qc/multiqc_fastqc.html",
        #"qc/multiqc_bam.html"


# 1 - trim reads
rule trimmomatic_pe:
    input:
        unpack(get_fastq)
    output:
        r1 = temp("trimmed/{sample}-{run}.1P.fastq.gz"),
        r2 = temp("trimmed/{sample}-{run}.2P.fastq.gz"),
        # reads where trimming entirely removed the mate
        r1_unpaired = temp("trimmed/{sample}-{run}.1U.fastq.gz"),
        r2_unpaired = temp("trimmed/{sample}-{run}.2U.fastq.gz")
    params:
        trimmer = ["LEADING:3","TRAILING:3","SLIDINGWINDOW:4:15","MINLEN:36","ILLUMINACLIP:{}:2:30:10".format(config["adapters"])]
    log:
        "logs/trimmomatic/{sample}-{run}.log"
    threads: 8
    wrapper:
        "0.73.0/bio/trimmomatic/pe"

# 2 - map reads
# !!! reads that became unpaired during trimming are not mapped!
rule map_reads:
    input:
        reads = ["trimmed/{sample}-{run}.1P.fastq.gz", "trimmed/{sample}-{run}.2P.fastq.gz"]
    output:
        temp("mapped/{sample}-{run}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}-{run}.log"
    params:
        index = config["ref_genome"],
        extra = r"-R '@RG\tID:{sample}-{run}\tSM:{sample}\tPL:ILLUMINA'",
        sort = "samtools",
        sort_order = "coordinate"
    threads: 8
    wrapper:
        "0.73.0/bio/bwa/mem"

# 3 - filter reads
# reads eliminated: read unmapped, mate unmapped, read fails platform/vendor check (-F 524)
rule samtools_view:
    input:
        "mapped/{sample}-{run}.sorted.bam"
    output:
        temp("filtered/{sample}-{run}.filtered.bam")
    params:
        "-b -F 524" # optional params string
    threads: 8
    wrapper:
        "0.73.0/bio/samtools/view"

# 4 - combine bam by sample
# combine all runs into one bam per sample
rule samtools_merge:
    input:
        get_sample_runs
    output:
        temp("merged/{sample}.merged.bam")
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@
    wrapper:
        "0.73.0/bio/samtools/merge"

# 5 - remove duplicates
# remove PCR and optical duplicate reads with Picard
rule mark_duplicates:
    input:
        "merged/{sample}.merged.bam"
    output:
        bam="dedup/{sample}.dedup.bam",
        metrics="dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=true",
        "CREATE_INDEX=true"
    wrapper:
        "0.73.0/bio/picard/markduplicates"


# QC rules need to be reviewed
# QCs
rule fastqc_raw:
    '''
    runs fastqc for both fastq files for each sample-run combo
    calling unpack(get_fastq) assigns each fastq to input.r1 & input.r2, respectively
    from input.fastq.gz fastqc creates input_fastqc.html & input_fastq.zip
    since the output files cannot be renamed during fastqc run, params infers their filenames
    and files are renamed in shell command
    '''
    input:
        unpack(get_fastq)
    output:
        r1_html = "qc/fastqc_raw/{sample}-{run}_R1_fastqc.html",
        r1_zip = "qc/fastqc_raw/{sample}-{run}_R1_fastqc.zip",
        r2_html = "qc/fastqc_raw/{sample}-{run}_R2_fastqc.html",
        r2_zip = "qc/fastqc_raw/{sample}-{run}_R2_fastqc.zip"
    params: # from input.fastq.gz fastqc creates input_fastqc.html & input_fastq.zip
        r1_html = lambda wildcards, input: input.r1.split('/')[-1].replace('.fastq.gz','_fastqc.html').replace('fq.gz', '_fastqc.html'),
        r1_zip = lambda wildcards, input: input.r1.split('/')[-1].replace('.fastq.gz', '_fastqc.zip').replace('fq.gz', '_fastqc.zip'),
        r2_html = lambda wildcards, input: input.r2.split('/')[-1].replace('.fastq.gz','_fastqc.html').replace('fq.gz', '_fastqc.html'),
        r2_zip = lambda wildcards, input: input.r2.split('/')[-1].replace('.fastq.gz', '_fastqc.zip').replace('fq.gz', '_fastqc.zip')
    threads: 8
    conda:
        "fastqc.yaml"
    log:
        "logs/fastqc_raw/{sample}-{run}.log"
    shell:
        '''
        fastqc -t {threads} -o qc/fastqc_raw {input.r1} {input.r2} 2> {log} \\
        && mv qc/fastqc_raw/{params.r1_html} {output.r1_html} \\
        && mv qc/fastqc_raw/{params.r1_zip} {output.r1_zip} \\
        && mv qc/fastqc_raw/{params.r2_html} {output.r2_html} \\
        && mv qc/fastqc_raw/{params.r2_zip} {output.r2_zip}
        '''

rule fastqc_trimmed:
    '''
    same as fastqc_raw but for trimmed fastq
    '''
    input:
        r1 = "trimmed/{sample}-{run}.1P.fastq.gz",
        r2 = "trimmed/{sample}-{run}.2P.fastq.gz"
    output:
        r1_html = "qc/fastqc_trimmed/{sample}-{run}.1P_fastqc.html",
        r1_zip = "qc/fastqc_trimmed/{sample}-{run}.1P_fastqc.zip",
        r2_html = "qc/fastqc_trimmed/{sample}-{run}.2P_fastqc.html",
        r2_zip = "qc/fastqc_trimmed/{sample}-{run}.2P_fastqc.zip"
    threads: 8
    conda:
        "fastqc.yaml"
    log:
        "logs/fastqc_trimmed/{sample}-{run}.log"
    shell:
        '''
        fastqc -t {threads} -o qc/fastqc_trimmed {input.r1} {input.r2} 2> {log}
        '''

rule multiqc_fastqc:
    '''
    runs multiqc on the fastqc reports for the raw and trimmed fastq in qc/fastqc_raw/ and qc/fastqc_trimmed/
    '''
    input:
        ["qc/fastqc_raw/{0}-{1}_R1_fastqc.html".format(c[0],c[1]) for c in SAMPLE_RUN_COMBS],
        ["qc/fastqc_trimmed/{0}-{1}.1P_fastqc.html".format(c[0],c[1]) for c in SAMPLE_RUN_COMBS],
    output:
        "qc/multiqc_fastqc.html"
    log:
        "logs/multiqc_fastqc.log"
    wrapper:
        "0.38.0/bio/multiqc"

# 5 - run QC on bam
rule samtools_stats_mapped:
    input:
        "mapped/{sample}-{run}.sorted.bam"
    output:
        "qc/mapped_stats/{sample}-{run}.txt"
    log:
        "logs/mapped_stats/{sample}-{run}.log"
    wrapper:
        "0.38.0/bio/samtools/stats"

rule samtools_stats_filtered:
    input:
        "filtered/{sample}-{run}.filtered.bam"
    output:
        "qc/filtered_stats/{sample}-{run}.filtered.txt"
    log:
        "logs/filtered_stats/{sample}-{run}.filtered.log"
    wrapper:
        "0.38.0/bio/samtools/stats"

rule multiqc_bam:
    input:
        ["qc/mapped_stats/{0}-{1}.txt".format(c[0],c[1]) for c in SAMPLE_RUN_COMBS],
        ["qc/filtered_stats/{0}-{1}.filtered.txt".format(c[0],c[1]) for c in SAMPLE_RUN_COMBS]
    output:
        "qc/multiqc_bam.html"
    log:
        "logs/multiqc_bam.log"
    wrapper:
        "0.38.0/bio/multiqc"

