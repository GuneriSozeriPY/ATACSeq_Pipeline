

configfile: "config/config.yaml"

include: "workflows/rules/common.smk"

rule all:
    input:
        #config["mysamples"],
        lambda w: allInput(config["mysamples"]),

include: "workflows/rules/unzip.smk"

include: "workflows/rules/fastqc.smk"