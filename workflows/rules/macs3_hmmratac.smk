rule macs3_cutoff_analysis:
    input:
        bam="results/{sample}/{sample}_{build}_sorted.blacklist-filtered.bam",
        idx="results/{sample}/{sample}_{build}_sorted.blacklist-filtered.bam.bai"
    output:
        cutoff="results/macs3_hmmratac/{sample}_{build}_hmmratac_cutoff_analysis.tsv"
    log:
        "logs/macs3_hmmratac/{sample}_{build}_analysis.log"
    params:
        outdir="results/macs3_hmmratac",
       name="{sample}_{build}_hmmratac",
        extra="-f BAMPE --cutoff-analysis-only" 
    conda:
        "../envs/macs3.yaml"
    shell:
        """
        (echo "Running Cutoff Analysis..." &&        
        # 1. Run MACS3
        macs3 hmmratac -i {input.bam} \
            --outdir {params.outdir} \
            -n {params.name} {params.extra})
        """

# --- Rule 2: Final HMMRATAC Calling (Logic embedded in Shell) ---
rule macs3_hmmratac_final:
    input:
        bam="results/{sample}/{sample}_{build}_sorted.blacklist-filtered.bam",
        idx="results/{sample}/{sample}_{build}_sorted.blacklist-filtered.bam.bai",
        cutoff="results/macs3_hmmratac/{sample}_{build}_hmmratac_cutoff_analysis.tsv"
    output:
        gapped_peak="results/macs3_hmmratac/{sample}_{build}_hmmratac_final_accessible_regions.gappedPeak",
    log:
        "logs/macs3_hmmratac/{sample}_{build}_calling.log"
    params:
        outdir="results/macs3_hmmratac/",
        name="{sample}_{build}_hmmratac_final",
        base_extra="-f BAMPE"
    conda:
        "../envs/macs3.yaml"
    shell:
        """
        (echo "Calculating dynamic cutoffs from {input.cutoff}..." &&
        
        # --- Embedded Python Script to Calculate Params ---
        # This ensures calculation happens AFTER the cutoff file exists
        DYNAMIC_PARAMS=$(python3 -c '
import sys
try:
    with open(sys.argv[1], "r") as f:
        lines = f.readlines()
    data = []
    for line in lines[1:]:
        p = line.strip().split()
        if len(p) >= 2: data.append((float(p[0]), int(p[1])))
    
    if not data:
        # Default fallback
        print("--upper 20 --lower 10 --prescan-cutoff 1.2")
    else:
        # Logic: Lower (-l) ~10k peaks
        l_val = int(min(data, key=lambda x: abs(x[1] - 10000))[0])

        # Logic: Upper (-u) < 1k peaks
        upper_cands = [x for x in data if x[1] < 1000]
        if upper_cands:
            u_val = int(max(upper_cands, key=lambda x: x[1])[0])
        else:
            u_val = int(data[0][0])

        # Logic: Prescan (-c) bottom
        c_val = data[-1][0]
        if c_val <= 1.0: c_val = 1.1  # Force > 1.0

        print(f"--upper {{u_val}} --lower {{l_val}} --prescan-cutoff {{c_val}}")
except:
    print("--upper 20 --lower 10 --prescan-cutoff 1.2")
' {input.cutoff}) &&

        echo "Dynamic Params Calculated: $DYNAMIC_PARAMS" &&
        
        # --- Run Final MACS3 ---
        macs3 hmmratac -i {input.bam} \
            --outdir {params.outdir} \
            -n {params.name} \
            {params.base_extra} \
            $DYNAMIC_PARAMS &&
            
        echo "Done" )
        """