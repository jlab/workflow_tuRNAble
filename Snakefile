import gzip
from scripts.utils import *


configfile: "settings.yaml"

localrules: split_rfam

rule download_reference_assembly:
    output:
        "{prefix}/Results/References/Assemblies/{species}/{assembly}.fa.gz"
    params:
        url=lambda wildcards: get_assembly_url(wildcards, config),
        proxy=config['system']['proxy']['https']
    log:
        "{prefix}/Logs/References/Assemblies/{species}/{assembly}.log"
    benchmark:
        "{prefix}/Benchmarks/References/Assemblies/{species}/{assembly}.benchmark"
    shell:
        "export https_proxy='{params.proxy}' 2> {log}"
        " && wget {params.url} -O - > {output} 2>> {log}"

rule download_rfam:
    output:
        "{prefix}/Results/References/Rfam/{version}/Rfam.cm.gz"
    params:
        url=lambda wildcards: get_rfam_url(wildcards, config),
        ftpproxy=config['system']['proxy']['https']
    log:
        "{prefix}/Logs/References/Rfam/{version}.log"
    benchmark:
        "{prefix}/Benchmarks/References/Rfam/{version}.benchmark"
    shell:
        "export https_proxy='{params.ftpproxy}' 2> {log}"
        " && wget {params.url} -O - > {output} 2>> {log}"

# we split the Rfam CM collection into individual families, as homology search
# in one go takes days AND did crash the process in my hands multiple times
rule split_rfam:
    output:
        directory("{prefix}/Results/References/Rfam/{version}/splitted/")
    input:
        "{prefix}/Results/References/Rfam/{version}/Rfam.cm.gz"
    log:
        "{prefix}/Logs/References/Rfam/split_rfam/{version}.log"
    benchmark:
        "{prefix}/Benchmarks/References/Rfam/split_rfam/{version}.benchmark"
    run:
        os.makedirs(output[0], exist_ok=True)
        with gzip.open(input[0]) as f:
            CM_started = False
            HMM_started = False
            accession = None

            single_model_lines = []
            for line in f.readlines():
                line = line.decode()
                single_model_lines.append(line)
                if line.startswith("ACC "):
                    accession = ' '.join(line.split()[1:]).strip()
                elif line.startswith('INFERNAL1/'):
                    CM_started = True
                elif line.startswith('HMMER3/'):
                    HMM_started = True
                elif line.strip() == "//":
                    if CM_started and HMM_started and (accession is not None):
                        fp = os.path.join(output[0], '%s.cm' % accession)
                        with open(fp, 'w') as w:
                            w.write(''.join(single_model_lines))
                        single_model_lines = []
                        CM_started = False
                        HMM_started = False

HUGE_MODELS = ['RF00017', 'RF01577', 'RF02543', 'RF01295', 'RF04128',
               'RF02540', 'RF02541',  # OOM with 8GB on BCF cluster
              ]
def get_mem(wildcards):
    if wildcards['family'] in HUGE_MODELS:
        return 120000  # in MB
    return 8000  # in MB
def get_partition(wildcards):
    if wildcards['family'] in HUGE_MODELS:
        return "idle"
    return "bcf"
def get_threads(wildcards):
    if wildcards['family'] in HUGE_MODELS:
        return 32
    return 4

rule find_ncRNAs_oneFamily:
    input:
        assembly="{prefix}/Results/References/Assemblies/{species}/{assembly}.fa.gz",
        cm="{prefix}/Results/References/Rfam/{rfam}/splitted/{family}.cm"
    output:
        tab="{prefix}/Results/ncRNAhits/{species}/{assembly}/splitted/{species}_assembly-{assembly}_rfam-{rfam}_{family}.tab",
        out="{prefix}/Results/ncRNAhits/{species}/{assembly}/splitted/{species}_assembly-{assembly}_rfam-{rfam}_{family}.out"
    log:
        "{prefix}/Logs/ncRNAhits/{species}/{assembly}/{rfam}/splitted/{family}.log"
    benchmark:
        "{prefix}/Benchmarks/ncRNAhits/{species}/{assembly}/{rfam}/splitted/{family}.benchmark"
    conda:
        "envs/all.yaml"
    resources:
        threads=get_threads,
        mem_mb=get_mem,
        partition=get_partition,
    shell:
        "cmsearch --cpu {resources.threads} --noali --tblout {output.tab} {input.cm} {input.assembly} > {output.out} 2> {log}"

def get_family_list(wildcards):
    fp_collection = "%s/Results/References/Rfam/%s/Rfam.cm.gz" % (wildcards.prefix, wildcards.rfam)
    families = []
    with gzip.open(fp_collection) as f:
        for line in f.readlines():
            line = line.decode()
            if line.startswith("ACC "):
                families.append(line.split()[1].strip())
    return ['%s/Results/ncRNAhits/%s/%s/splitted/%s_assembly-%s_rfam-%s_%s.tab' % (
        wildcards.prefix, wildcards.species, wildcards.assembly,
        wildcards.species, wildcards.assembly, wildcards.rfam, f) for f in sorted(list(set(families)))]
rule find_ncRNAs:
    input:
        assembly="{prefix}/Results/References/Assemblies/{species}/{assembly}.fa.gz",
        family_results=lambda wildcards: get_family_list(wildcards)#  lambda wildcards: expand("{prefix}/Results/References/Rfam/{rfam}/splitted/{family}.tab", family=get_family_list("%s/Results/References/Rfam/%s/Rfam.cm.gz" % (wildcards.prefix, wildcards.rfam))),
    output:
        tab="{prefix}/Results/ncRNAhits/{species}/{assembly}/{species}_assembly-{assembly}_rfam-{rfam}.tab",
        out="{prefix}/Results/ncRNAhits/{species}/{assembly}/{species}_assembly-{assembly}_rfam-{rfam}.out"
    log:
        "{prefix}/Logs/ncRNAhits/{species}/{assembly}/{rfam}.log"
    benchmark:
        "{prefix}/Benchmarks/ncRNAhits/{species}/{assembly}/{rfam}.benchmark"
    conda:
        "envs/all.yaml"
    shell:
        "cat {input} > {output} 2> {log}"
