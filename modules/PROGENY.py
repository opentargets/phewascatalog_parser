import pandas as pd
import argparse
import gzip
import logging
import json
from pyspark import SparkContext, SparkFiles
from pyspark.sql import SparkSession
from pyspark.sql.functions import *
from pyspark.sql.types import *

# Pathway -> Perturbed Targets
# https://drive.google.com/drive/folders/1L5Y_umEZiccWJnXiiaYMNKUKYTjnp3ZU
PATHWAY_TARGET_MAP = {
    'Androgen': ['AR'],
    'EGFR': ['EGFR'],
    'Estrogen': ['ESR1'],
    'Hypoxia': ['HIF1A'],
    'JAK.STAT': ['JAK1', 'JAK2', 'STAT1', 'STAT2'],
    'MAPK': ['MAPK2K1', 'MAP2K2', 'RAF1'],
    'NFkB': ['TLR4', 'NKFB1', 'RELA'],
    'PI3K': ['PIK3CA', 'PI3K(Class1)'],
    'TGFb': ['TGFBR1', 'TGFBR2'],
    'TNFa': ['TNFRSF1A'],
    'Trail': ['TNFSF10', 'BCL2', 'BCL-XL', 'BCL-W', 'MCL1'],
    'VEGF': ['VEGFR'],
    'WNT': ['WNT3A', 'GSK3A', 'GSK3B'],
    'p53': ['TP53']
}

# === Pathway -> Reactome Pathway ID ===
#TODO: 1. For Hypoxia and Trail pathway mapping to be established and description to be updated
#TODO: 2. Androgen, Estrogen and WNT pathways should be referenced with https://www.biorxiv.org/content/10.1101/532739v1
#TODO: 3. MAPK and p53 are mapped to the same pathway ID
PATHWAY_REACTOME_MAP = {
    'Androgen': 'R-HSA-8940973:RUNX2 regulates osteoblast differentiation',
    'EGFR': 'R-HSA-8856828:Clathrin-mediated endocytosis',
    'Estrogen': 'R-HSA-8939902:Regulation of RUNX2 expression and activity',
    'Hypoxia': 'R-HSA-123456:XXX Pathway desc to be updated',
    'JAK.STAT': 'R-HSA-6785807:Interleukin-4 and 13 signaling',
    'MAPK': 'R-HSA-2559580:Oxidative Stress Induced Senescence',
    'NFkB': 'R-HSA-9020702:Interleukin-1 signaling',
    'PI3K': 'R-HSA-8853659:RET signaling',
    'TGFb': 'R-HSA-2173788:Downregulation of TGF-beta receptor signaling',
    'TNFa': 'R-HSA-5357956:TNFR1-induced NFkappaB signaling pathway',
    'Trail': 'R-HSA-123456:XXX Pathway desc to be updated',
    'VEGF': 'R-HSA-1234158:Regulation of gene expression by Hypoxia-inducible Factor',
    'WNT': 'R-HSA-381340:Transcriptional regulation of white adipocyte differentiation',
    'p53': 'R-HSA-2559580:Oxidative Stress Induced Senescence'
}

# === These symbols are secondary/generic/typo that need updating ===
PROGENY_SYMBOL_MAPPING = {
    'NKFB1': 'NFKB1',
    'MAPK2K1': 'PRKMK1',
    'PI3K(Class1)': 'PIK3CA',
    'VEGFR': 'KDR',
    'BCL-W': 'BCL2L2',
    'BCL-XL': 'BCL2L1'
}

class progenyEvidenceGenerator():
    def __init__(self, inputFile, mappingStep):
        # Create spark session     
        self.spark = SparkSession.builder \
                .appName('evidence_builder') \
                .getOrCreate()

        # Initialize mapping variables
        self.mappingStep = mappingStep
    
        # Initialize input files
        self.inputFile = inputFile
        self.dataframe = None

    def writeEvidenceFromSource(self):
        '''
        Processing of the input file to build all the evidences from its data
        Returns:
            evidences (array): Object with all the generated evidences strings from source file
        '''
        # Read input file
        self.dataframe = self.spark.read \
                                .option("header", "true") \
                                .option("delimiter", "\t") \
                                .option("inferSchema", "true") \
                                .csv(self.inputFile)

        # Mapping step
        if self.mappingStep:
            self.dataframe = self.cancer2EFO()
        
        self.dataframe = self.pathway2Reactome()

        # Build evidence strings per row
        evidences = self.dataframe.rdd \
            .map(progenyEvidenceGenerator.parseEvidenceString) \
            .collect() # list of dictionaries
        
        return evidences
    
    def cancer2EFO(self):
        diseaseMappingsFile = self.spark \
                        .read.csv("resources/cancer2EFO_mappings.tsv", sep=r'\t', header=True) \
                        .select("Cancer_type_acronym", "EFO_id") \
                        .withColumnRenamed("Cancer_type_acronym", "Cancer_type")

        self.dataframe = self.dataframe.join(
            diseaseMappingsFile,
            on="Cancer_type",
            how="inner"
        )

        return self.dataframe
    
    def pathway2Reactome(self):
        pathwayMappingsFile = self.spark \
                        .read.csv("resources/pathway2Reactome_mappings.tsv", sep=r'\t', header=True) \
                        .withColumnRenamed("pathway", "Pathway")
        
        self.dataframe = self.dataframe \
                .join(
                    pathwayMappingsFile,
                    on="Pathway",
                    how="inner"
                )
        self.dataframe = self.dataframe \
                .withColumn("target", explode("target"))
        #print(self.dataframe.first())

        return self.dataframe

    @staticmethod
    def parseEvidenceString(row):
        try:
            evidence = {
                "datasourceId" : "progeny",
                "datatypeId" : "affected_pathway",
                "diseaseFromSource" : row["Cancer_type"],
                "diseaseFromSourceMappedId" : row["EFO_id"],
                "resourceScore" : row["P.Value"],
                "pathwayName" : row["description"],
                "pathwayId" : row["reactomeId"],
                "targetFromSourceId" : row["target"] # TO-DO: mapping corrections with PROGENY_SYMBOL_MAPPING
            }
            return evidence
        except Exception as e:
            raise        

def main():
    # Initiating parser
    parser = argparse.ArgumentParser(description=
    "This script generates evidences for the PROGENy data source.")

    parser.add_argument("-i", "--inputFile", required=True, type=str, help="Input .csv file with the table containing association details.")
    parser.add_argument("-o", "--outputFile", required=True, type=str, help="Name of the json output file containing the evidence strings.")
    parser.add_argument("-m", "--mappingStep", required=False, type=bool, default=True, help="State whether to run the disease to EFO term mapping step or not.")

    # Parsing parameters
    args = parser.parse_args()

    inputFile = args.inputFile
    outputFile = args.outputFile
    mappingStep = args.mappingStep

    # Initialize logging:
    logging.basicConfig(
    filename='evidence_builder.log',
    level=logging.INFO,
    format='%(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    )

    # Initialize evidence builder object
    evidenceBuilder = progenyEvidenceGenerator(inputFile, mappingStep)

    # Writing evidence strings into a json file
    evidences = evidenceBuilder.writeEvidenceFromSource()

    with gzip.open(outputFile, "wt") as f: # TO-DO: export in .gz
        for evidence in evidences:
            json.dump(evidence, f)
            f.write('\n')
    logging.info(f"Evidence strings saved into {outputFile}. Exiting.")

if __name__ == '__main__':
    main()