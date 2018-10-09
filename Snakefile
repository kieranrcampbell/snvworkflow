configfile: "config.yaml"

ac_names = ["covmat.mtx", "refmat.mtx", "altmat.mtx", "-barcode_index.csv", "-region_index.csv"]

ac_files = expand(config['allelecount_dir'] + "/" + config['id'] + "{file}",file=ac_names)

rule all:
    input:
        config['cnv_clone_mat'],
        config['cov'],
        config['ref']

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
        

rule snvworkflow:
    input:
        clone_cnv = config['input_clone_cnv'],
        acf = ac_files
    output:
        cnv_clone_mat = config['cnv_clone_mat'],
        cov = config['cov'],
        ref = config['ref']
    shell:
        "Rscript scripts/snvworkflow.R \
        --input_clone_cnv {input.clone_cnv} \
        --input_cov {ac_files[0]} \
        --input_ref {ac_files[1]} \
        --input_alt {ac_files[2]} \
        --input_barcode_index {ac_files[3]} \
        --input_region_index {ac_files[4]} \
        --output_cnv_clone {output.cnv_clone_mat} \
        --output_cov {output.cov} \
        --output_ref {output.ref}"
