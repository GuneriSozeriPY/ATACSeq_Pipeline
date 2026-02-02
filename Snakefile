configfile: "config/config.yaml"

include: "workflows/rules/common.smk"

wildcard_constraints:
    build=config["genome"]["build"],

rule all:
    input:
        lambda w: allInput(config["genome"]["build"], config["meta"], config.get("conditions"))

include: "workflows/rules/prepare_genome.smk"

include: "workflows/rules/fastqc.smk"

include: "workflows/rules/fastp.smk"

include: "workflows/rules/alignment.smk"

include: "workflows/rules/macs2_idr.smk"


