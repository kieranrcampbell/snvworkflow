configfile: "config.yaml"

ac_names = ["covmat.mtx", "refmat.mtx", "altmat.mtx", "-barcode_index.csv", "-region_index.csv"]

ac_files = expand(config['allelecount_dir'] + "/" + config['id'] + "{file}",file=ac_names)

rule all:
    input:
        ac_files

rule run_allelecount:
    params:
        output = config['allelecount_dir'] + "/" + config['id'],
	snps = config['snps'],
	barcodes = config['barcodes'],
	bam10X = config['bam10X']
    input:
        config['bam10X'], config['snps'], config['barcodes']
    output:
        ac_files
    shell:
         "python scripts/scAlleleCount.py -v --snps {params.snps} --barcodes {params.barcodes} --output-prefix {params.output} {params.bam10X}"
        

