configfile: "config/config.yaml"

include: "workflows/rules/common.smk"

wildcard_constraints:
    build=config["genome"]["build"],

rule all:
    input:
        lambda w: allInput(config["genome"]["build"], config["meta"]),

include: "workflows/rules/prepare_genome.smk"

include: "workflows/rules/fastqc.smk"

include: "workflows/rules/fastp.smk"

include: "workflows/rules/alignment.smk"

include: "workflows/rules/macs2.smk"