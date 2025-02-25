

def allInput(metadata):
    inputlist = []
    for sample in metadata:
        inputlist.append(f"results/{sample}_fastqc.html")
    return inputlist