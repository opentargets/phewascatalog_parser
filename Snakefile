from datetime import datetime
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

GS = GSRemoteProvider()
HTTP = HTTPRemoteProvider()

# --- Settings --- #
# Current date in YYYY-MM-DD format:
timeStamp = datetime.now().strftime("%Y-%m-%d")

# Configuration is read from the config yaml:
configfile: 'configuration.yaml'
logFile = f"{config['global']['logDir']}/evidence_parser-{timeStamp}.log"

# --- All rules --- #
rule all:
    input:
        GS.remote(f"{config['ClinGen']['outputBucket']}/ClinGen-{timeStamp}.json.gz"),
        GS.remote(f"{config['PheWAS']['outputBucket']}/phewas_catalog-{timeStamp}.json.gz"),
        GS.remote(f"{config['SLAPEnrich']['outputBucket']}/slapenrich-{timeStamp}.json.gz"),
        GS.remote(f"{config['Gene2Phenotype']['outputBucket']}/gene2phenotype-{timeStamp}.json.gz"),
        GS.remote(f"{config['CRISPR']['outputBucket']}/crispr-{timeStamp}.json.gz"),
        GS.remote(f"{config['SysBio']['outputBucket']}/sysbio-{timeStamp}.json.gz"),
        directory(GS.remote(f"{config['EPMC']['outputBucket']}/epmc-{timeStamp}")),
        directory(GS.remote(gs://genetics-portal-analysis/l2g-platform-export/data/genetics_portal_evidence.json.gz)),
        GS.remote(f"{config['PROGENy']['outputBucket']}/progeny-{timeStamp}.json.gz"),
        GS.remote(f"{config['intOGen']['outputBucket']}/intogen-{timeStamp}.json.gz"),
        GS.remote(f"{config['PanelApp']['outputBucket']}/genomics_england-{timeStamp}.json.gz"),
        GS.remote(f"{config['Phenodigm']['outputBucket']}/phenodigm-{timeStamp}.json.gz"),
        GS.remote(f"{config['Orphanet']['outputBucket']}/Orphanet-{timeStamp}")

# --- Auxiliary Rules --- #
## help                     : prints help comments for Snakefile
rule help:
    input: "Snakefile"
    shell:
        "sed -n 's/^##//p' {input}"

## clean                    : prints help comments for Snakefile
rule clean:
    shell:
        "rm -rf tmp"

# --- Data sources parsers --- #
## clingen                  : processes the Gene Validity Curations table from ClinGen
rule clingen:
    input:
        f"tmp/ClinGen-Gene-Disease-Summary-{timeStamp}.csv"
    output:
        evidenceFile=GS.remote(f"{config['ClinGen']['outputBucket']}/ClinGen-{timeStamp}.json.gz"),
        unmappedDiseases=f"tmp/unmappedDiseases/clingen_unmapped_diseases-{timeStamp}.lst"
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/ClinGen.py \
        --input_file {input} \
        --output_file {output.evidenceFile}
        """

## geneticsPortal           : processes lead variants from the Open Targets Genetics portal on a Dataproc cluster
rule geneticsPortal:
    shell:
        """
        gcloud dataproc clusters create \
            snakemake-cluster-l2g-data \
            --image-version=1.4 \
            --properties=spark:spark.debug.maxToStringFields=100,spark:spark.executor.cores=31,spark:spark.executor.instances=1 \
            --master-machine-type=n1-standard-32 \
            --master-boot-disk-size=1TB \
            --zone=europe-west1-d \
            --single-node \
            --max-idle=5m \
            --region=europe-west1 \
            --project=open-targets-genetics
        
        gcloud dataproc jobs submit pyspark \
            --cluster=snakemake-cluster-l2g-data \
            --project=open-targets-genetics \
            --region=europe-west1 \
            modules/GeneticsPortal.py -- \
            --locus2gene gs://genetics-portal-data/l2g/200127 \
            --toploci gs://genetics-portal-data/v2d/200207/toploci.parquet \
            --study gs://genetics-portal-data/v2d/200207/studies.parquet \
            --threshold 0.05 \
            --variantIndex gs://genetics-portal-data/variant-annotation/190129/variant-annotation.parquet  \
            --ecoCodes gs://genetics-portal-data/lut/vep_consequences.tsv \
            --outputFile gs://genetics-portal-analysis/l2g-platform-export/data/genetics_portal_evidence.json.gz
        """

## phewas                   : processes phenome-wide association studies data from PheWAS
rule phewas:
    input:
        inputFile=f"tmp/phewas_catalog-{timeStamp}.csv",
        consequencesFile=f"tmp/phewas_w_consequences-{timeStamp}.csv",
        diseaseMapping=f"tmp/phewascat_mappings-{timeStamp}.tsv",
    params:
        genesSet=f"{config['global']['genesHGNC']}"
    output:
        evidenceFile=GS.remote(f"{config['PheWAS']['outputBucket']}/phewas_catalog-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        '''
        python modules/PheWAS.py \
            --inputFile {input.inputFile} \
            --consequencesFile {input.consequencesFile} \
            --diseaseMapping {input.diseaseMapping} \
            --genesSet {params.genesSet} \
            --outputFile {output.evidenceFile}
        '''

## slapenrich               : processes cancer-target evidence strings derived from SLAPenrich
rule slapenrich:
    input:
        inputFile=f"tmp/slapenrich-{timeStamp}.csv",
        diseaseMapping=f"tmp/cancer2EFO_mappings-{timeStamp}.tsv"
    output:
        evidenceFile=GS.remote(f"{config['SLAPEnrich']['outputBucket']}/slapenrich-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/SLAPEnrich.py \
            --inputFile {input.inputFile} \
            --diseaseMapping {input.diseaseMapping} \
            --outputFile {output.evidenceFile} \
        """

## gene2Phenotype           : processes four gene panels from Gene2Phenotype
rule gene2Phenotype:
    input:
        ddPanel=f"tmp/DDG2P-{timeStamp}.csv.gz",
        eyePanel=f"tmp/EyeG2P-{timeStamp}.csv.gz",
        skinPanel=f"tmp/SkinG2P-{timeStamp}.csv.gz",
        cancerPanel=f"tmp/CancerG2P-{timeStamp}.csv.gz"
    output:
        evidenceFile=GS.remote(f"{config['Gene2Phenotype']['outputBucket']}/gene2phenotype-{timeStamp}.json.gz"),
        unmappedDiseases=f"tmp/unmappedDiseases/gene2phenotype_unmapped_diseases-{timeStamp}.txt"
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/Gene2Phenotype.py \
            --dd_panel {input.ddPanel} \
            --eye_panel {input.eyePanel} \
            --skin_panel {input.skinPanel} \
            --cancer_panel {input.cancerPanel} \
            --output_file {output.evidenceFile}
        """

## crispr                   : processes cancer therapeutic targets using CRISPR–Cas9 screens
rule crispr:
    input:
        evidenceFile=f"tmp/crispr_evidence-{timeStamp}.csv",
        descriptionsFile=f"tmp/crispr_descriptions-{timeStamp}.tsv",
        cellTypesFile=f"tmp/crispr_cell_lines-{timeStamp}.tsv"
    output:
        evidenceFile=GS.remote(f"{config['CRISPR']['outputBucket']}/crispr-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/CRISPR.py \
            --evidence_file {input.evidenceFile} \
            --descriptions_file {input.descriptionsFile} \
            --cell_types_file {input.cellTypesFile} \
            --output_file {output.evidenceFile}
        """

## progeny                  : processes gene expression data from TCGA derived from PROGENy
rule progeny:
    input:
        inputFile=f"tmp/progeny_normalVStumor_opentargets-{timeStamp}.txt",
        diseaseMapping=f"tmp/progeny_cancer2EFO_mappings-{timeStamp}.tsv",
        pathwayMapping=f"tmp/pathway2Reactome_mappings-{timeStamp}.tsv"
    output:
        evidenceFile=GS.remote(f"{config['PROGENy']['outputBucket']}/progeny-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/PROGENY.py \
            --inputFile {input.inputFile} \
            --diseaseMapping {input.diseaseMapping} \
            --pathwayMapping {input.pathwayMapping} \
            --outputFile {output.evidenceFile}
        """

## phenodigm                : processes target-disease evidence querying the IMPC SOLR API
rule phenodigm:
    output:
        evidenceFile=GS.remote(f"{config['Phenodigm']['outputBucket']}/phenodigm-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/PhenoDigm.py \
            --cache-dir phenodigm_cache \
            --output {output.evidenceFile}
        """

## sysbio                   : processes key driver genes for specific diseases that have been curated from Systems Biology papers
rule sysbio:
    input:
        evidenceFile=f"tmp/sysbio_evidence-{timeStamp}.tsv",
        studyFile=f"tmp/sysbio_publication_info-{timeStamp}.tsv"
    output:
        evidenceFile=GS.remote(f"{config['SysBio']['outputBucket']}/sysbio-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/SystemsBiology.py \
            --evidenceFile {input.evidenceFile} \
            --studyFile {input.studyFile} \
            --outputFile {output.evidenceFile}
        """

## panelApp                 : processes gene panels data curated by Genomics England
rule panelApp:
    input:
        inputFile=f"tmp/panelapp_gene_panels-{timeStamp}.tsv"
    output:
        evidenceFile = GS.remote(f"{config['PanelApp']['outputBucket']}/genomics_england-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/PanelApp.py \
            --inputFile {input.inputFile} \
            --outputFile {output.evidenceFile}
        """

## intogen                  : processes cohorts and driver genes data from intOGen
rule intogen:
    input:
        inputGenes = f"tmp/Compendium_Cancer_Genes-{timeStamp}.tsv",
        inputCohorts = f"tmp/cohorts-{timeStamp}.tsv",
        diseaseMapping=f"tmp/intogen_cancer2EFO_mappings-{timeStamp}.tsv"
    output:
        evidenceFile = GS.remote(f"{config['intOGen']['outputBucket']}/intogen-{timeStamp}.json.gz")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/IntOGen.py \
            --inputGenes {input.inputGenes} \
            --inputCohorts {input.inputCohorts} \
            --diseaseMapping {input.diseaseMapping} \
            --outputFile {output.evidenceFile}
        """

## epmc                     : processes target/disease evidence strings from ePMC cooccurrence files
rule epmc:
    input:
        inputCooccurences = directory(f"tmp/epmc_cooccurrences-{timeStamp}")
    output:
        evidenceFile = directory(GS.remote(f"{config['EPMC']['outputBucket']}/epmc-{timeStamp}"))
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/EPMC.py \
            --cooccurrenceFile {input.inputCooccurences} \
            --output {output.evidenceFile} \
            --local
        """


## Orphanet                 : Processing disease/target evidence from Orphanet
rule orphanet:
    input:
        HTTP.remote(config['Orphanet']['webSource'])
    output:
        GS.remote(f"{config['Orphanet']['outputBucket']}/Orphanet-{timeStamp}")
    log:
        GS.remote(logFile)
    shell:
        """
        python modules/Orphanet.py \
            --input_file {input} \
            --output_file {output} \
            --local
        """

# --- Fetching input data and uploading to GS --- #
## fetchClingen             : fetches the Gene Validity Curations table from ClinGen
rule fetchClingen:
    params:
        webSource=config['ClinGen']['webSource']
    output:
        bucket=GS.remote(f"{config['ClinGen']['inputBucket']}/ClinGen-Gene-Disease-Summary-{timeStamp}.csv"),
        local=f"tmp/ClinGen-Gene-Disease-Summary-{timeStamp}.csv"
    log:
        GS.remote(logFile)
    shell:
        """
        curl {params.webSource} > {output.local}
        gsutil cp {output.local} {output.bucket}
        """

## fetchPhewas              : fetches the PheWAS data and an enriched table from GS and a disease mapping look-up
rule fetchPhewas:
    input:
        inputFile=GS.remote(f"{config['PheWAS']['inputBucket']}/phewas-catalog-19-10-2018.csv"),
        consequencesFile=GS.remote(f"{config['PheWAS']['inputBucket']}/phewas_w_consequences.csv")
    params:
        diseaseMapping=config['PheWAS']['diseaseMapping'],
    output:
        inputFile=f"tmp/phewas_catalog-{timeStamp}.csv",
        consequencesFile=f"tmp/phewas_w_consequences-{timeStamp}.csv",
        diseaseMapping=f"tmp/phewascat_mappings-{timeStamp}.tsv"
    log:
        GS.remote(logFile)
    shell:
        """
        curl {params.diseaseMapping} > {output.diseaseMapping}
        gsutil cp {input.inputFile} {output.inputFile}
        gsutil cp {input.consequencesFile} {output.consequencesFile}
        """

## fetchSlapenrich          : fetches SLAPenrich table from GS
rule fetchSlapenrich:
    input:
        inputFile=GS.remote(f"{config['SLAPEnrich']['inputBucket']}/slapenrich_opentargets-21-12-2017.tsv")
    params:
        diseaseMapping=config['SLAPEnrich']['diseaseMapping']
    output:
        inputFile=f"tmp/slapenrich-{timeStamp}.csv",
        diseaseMapping=f"tmp/cancer2EFO_mappings-{timeStamp}.tsv"
    log:
        GS.remote(logFile)
    shell:
        """
        cp {params.diseaseMapping} {output.diseaseMapping}
        gsutil cp {input.inputFile} {output.inputFile}
        """

## fetchGene2Phenotype      : fetches four gene panels downloaded from Gene2Phenotype
rule fetchGene2Phenotype:
    params:
        webSource_dd_panel=config['Gene2Phenotype']['webSource_dd_panel'],
        webSource_eye_panel=config['Gene2Phenotype']['webSource_eye_panel'],
        webSource_skin_panel=config['Gene2Phenotype']['webSource_skin_panel'],
        webSource_cancer_panel=config['Gene2Phenotype']['webSource_cancer_panel']

    output:
        ddBucket=GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/DDG2P-{timeStamp}.csv.gz"),
        ddLocal=f"tmp/DDG2P-{timeStamp}.csv.gz",
        eyeBucket=GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/EyeG2P-{timeStamp}.csv.gz"),
        eyeLocal=f"tmp/EyeG2P-{timeStamp}.csv.gz",
        skinBucket=GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/SkinG2P-{timeStamp}.csv.gz"),
        skinLocal=f"tmp/SkinG2P-{timeStamp}.csv.gz",
        cancerBucket=GS.remote(f"{config['Gene2Phenotype']['inputBucket']}/CancerG2P-{timeStamp}.csv.gz"),
        cancerLocal=f"tmp/CancerG2P-{timeStamp}.csv.gz"
    log:
        GS.remote(logFile)
    shell:
        """
        curl {params.webSource_dd_panel} > {output.ddLocal}
        gsutil cp {output.ddLocal} {output.ddBucket}
        curl {params.webSource_eye_panel} > {output.eyeLocal}
        gsutil cp {output.eyeLocal} {output.eyeBucket}
        curl {params.webSource_skin_panel} > {output.skinLocal}
        gsutil cp {output.skinLocal} {output.skinBucket}
        curl {params.webSource_cancer_panel} > {output.cancerLocal}
        gsutil cp {output.cancerLocal} {output.cancerBucket}
        """

## fetchCrispr              : fetches three tables from GS
rule fetchCrispr:
    input:
        evidenceFile=GS.remote(f"{config['CRISPR']['inputBucket']}/crispr_evidence.tsv"),
        descriptionsFile=GS.remote(f"{config['CRISPR']['inputBucket']}/crispr_descriptions.tsv"),
        cellTypesFile=GS.remote(f"{config['CRISPR']['inputBucket']}/crispr_cell_lines.tsv")
    output:
        evidenceFile=f"tmp/crispr_evidence-{timeStamp}.csv",
        descriptionsFile=f"tmp/crispr_descriptions-{timeStamp}.tsv",
        cellTypesFile=f"tmp/crispr_cell_lines-{timeStamp}.tsv"
    log:
        GS.remote(logFile)
    shell:
        """
        gsutil cp {input.evidenceFile} {output.evidenceFile}
        gsutil cp {input.descriptionsFile} {output.descriptionsFile}
        gsutil cp {input.cellTypesFile} {output.cellTypesFile}
        """

## fetchProgeny             : fetches PROGENy table from GS, and two disease-mapping and pathway-mapping look-up tables
rule fetchProgeny:
    input:
        inputFile=GS.remote(f"{config['PROGENy']['inputBucket']}/progeny_normalVStumor_opentargets.txt"),
        diseaseMapping=config['PROGENy']['diseaseMapping'],
        pathwayMapping=config['PROGENy']['pathwayMapping']
    output:
        inputFile=f"tmp/progeny_normalVStumor_opentargets-{timeStamp}.txt",
        diseaseMapping=f"tmp/progeny_cancer2EFO_mappings-{timeStamp}.tsv", # solve ambiguity between progeny and slapenrich
        pathwayMapping=f"tmp/pathway2Reactome_mappings-{timeStamp}.tsv"
    log:
        GS.remote(logFile)
    shell:
        """
        cp {input.diseaseMapping} {output.diseaseMapping}
        cp {input.pathwayMapping} {output.pathwayMapping}
        gsutil cp {input.inputFile} {output.inputFile}
        """

## fetchSysbio              : fetches evidence data and study level information from GS
rule fetchSysbio:
    input:
        evidenceFile=GS.remote(f"{config['SysBio']['inputBucket']}/sysbio_evidence-31-01-2019.tsv"),
        studyFile=GS.remote(f"{config['SysBio']['inputBucket']}/sysbio_publication_info_nov2018.tsv")
    output:
        evidenceFile=f"tmp/sysbio_evidence-{timeStamp}.tsv",
        studyFile=f"tmp/sysbio_publication_info-{timeStamp}.tsv"
    log:
        GS.remote(logFile)
    shell:
        """
        gsutil cp {input.evidenceFile} {output.evidenceFile}
        gsutil cp {input.studyFile} {output.studyFile}
        """

## fetchPanelApp            : fetches gene panels data from GS
rule fetchPanelApp:
    input:
        inputFile = GS.remote(f"{config['PanelApp']['inputBucket']}/All_genes_20200928-1959.tsv")
    output:
        inputFile=f"tmp/panelapp_gene_panels-{timeStamp}.tsv"
    log:
        GS.remote(logFile)
    shell:
        """
        gsutil cp {input.inputFile} {output.inputFile}
        """

## fetchIntogen             : fetches cohorts and driver genes from GS
rule fetchIntogen:
    input:
        inputGenes = GS.remote(f"{config['intOGen']['inputBucket']}/Compendium_Cancer_Genes.tsv"),
        inputCohorts = GS.remote(f"{config['intOGen']['inputBucket']}/cohorts.tsv"),
        diseaseMapping=config['intOGen']['diseaseMapping']
    output:
        inputGenes = f"tmp/Compendium_Cancer_Genes-{timeStamp}.tsv",
        inputCohorts = f"tmp/cohorts-{timeStamp}.tsv",
        diseaseMapping=f"tmp/intogen_cancer2EFO_mappings-{timeStamp}.tsv"
    log:
        GS.remote(logFile)
    shell:
        """
        gsutil cp {input.inputGenes} {output.inputGenes}
        gsutil cp {input.inputCohorts} {output.inputCohorts}
        gsutil cp {input.diseaseMapping} {output.diseaseMapping}
        """

## fetchEpmc                : fetches the partioned parquet files with the ePMC cooccurrences
rule fetchEpmc:
    input:
        inputCooccurences = GS.remote(config['EPMC']['inputBucket'])
    output:
        inputCooccurences = directory(f"tmp/epmc_cooccurrences-{timeStamp}")
    log:
        GS.remote(logFile)
    shell:
        """
        gsutil cp -r {input.inputCooccurences} {output.inputCooccurences}
        """
