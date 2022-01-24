# only variable needed to change


PROJECT_PATH='/Users/yiquan/PycharmProjects/HA_egg_passage'
FQ = PROJECT_PATH + '/result/{SAMPLENAME}.assembled.fastq'


SAMPLENAMES, = glob_wildcards(FQ)

REF=PROJECT_PATH + '/ref/WT_seq.fasta'

RESULT_PATH = PROJECT_PATH + '/result/{SAMPLENAME}'

TABLE = RESULT_PATH + '_count.tsv'

rule all:
    input:
        expand(TABLE, SAMPLENAME=SAMPLENAMES)


rule fq2count:
    input:
        FQ
    params:
        REF_FA=REF
    output:
        TABLE
    shell:
        'python fastq2count.py {input} {params} {output}'
