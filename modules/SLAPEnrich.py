import argparse
import gzip
import logging
import json
from pyspark import SparkContext
from pyspark.sql import SparkSession
from pyspark.sql.functions import *
from pyspark.sql.types import *

class SLAPEnrichEvidenceGenerator():
    def __init__(self):
        # Create spark session     
        self.spark = SparkSession.builder \
                .appName('SLAPEnrich') \
                .getOrCreate()
    
        # Initialize source table
        self.dataframe = None

    def generateEvidenceFromSource(self, inputFile, skipMapping):
        '''
        Processing of the input file to build all the evidences from its data
        Returns:
            evidences (array): Object with all the generated evidences strings from source file
        '''
        # Read input file
        self.dataframe = self.spark \
                        .read.csv(inputFile, sep=r'\t', header=True, inferSchema=True) \
                        .select("ctype", "gene", "pathway", "SLAPEnrichPval") \
                        .withColumnRenamed("ctype", "Cancer_type_acronym") \
                        .withColumnRenamed("SLAPEnrichPval", "pval") \
                        .withColumn("pathwayId", split(col("pathway"), ": ").getItem(0)) \
                        .withColumn("pathwayDescription", split(col("pathway"), ": ").getItem(1))

        # Filter by p-value
        self.dataframe = self.dataframe.filter(col("pval") < 1e-4) 

        # Mapping step
        if not skipMapping:
            self.dataframe = self.cancer2EFO()

        # Build evidence strings per row
        evidences = self.dataframe.rdd \
            .map(SLAPEnrichEvidenceGenerator.parseEvidenceString) \
            .collect() # list of dictionaries
        
        return evidences
    
    def cancer2EFO(self):
        diseaseMappingsFile = self.spark \
                        .read.csv("resources/cancer2EFO_mappings.tsv", sep=r'\t', header=True) \
                        .select("Cancer_type_acronym", "EFO_id") \

        self.dataframe = self.dataframe.join(
            diseaseMappingsFile,
            on="Cancer_type_acronym",
            how="inner"
        )

        return self.dataframe

    @staticmethod
    def parseEvidenceString(row):
        try:
            evidence = {
                "datasourceId" : "slapenrich",
                "datatypeId" : "affected_pathway",
                "resourceScore" : row["pval"],
                "pathwayName" : row["pathwayDescription"],
                "pathwayId" : row["pathwayId"],
                "targetFromSourceId" : row["gene"],
                "diseaseFromSource" : row["Cancer_type_acronym"]
            }
            if "EFO_id" in row:
                evidence["diseaseFromSourceMappedId"] = row["EFO_id"]
            return evidence
        except Exception as e:
            raise        

def main():
    # Initiating parser
    parser = argparse.ArgumentParser(description=
    "This script generates evidences for the SLAPEnrich data source.")

    parser.add_argument("-i", "--inputFile", required=True, type=str, help="Input source .tsv file.")
    parser.add_argument("-o", "--outputFile", required=True, type=str, help="Name of the evidence compressed JSON file containing the evidence strings.")
    parser.add_argument("-s", "--skipMapping", required=False, action="store_true", help="State whether to skip the disease to EFO mapping step.")

    # Parsing parameters
    args = parser.parse_args()

    inputFile = args.inputFile
    outputFile = args.outputFile
    skipMapping = args.skipMapping

    # Initialize logging:
    logging.basicConfig(
    filename='evidence_builder.log',
    level=logging.INFO,
    format='%(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    )

    # Initialize evidence builder object
    evidenceBuilder = SLAPEnrichEvidenceGenerator()

    # Writing evidence strings into a json file
    evidences = evidenceBuilder.generateEvidenceFromSource(inputFile, skipMapping)

    with gzip.open(outputFile, "wt") as f:
        for evidence in evidences:
            json.dump(evidence, f)
            f.write('\n')
    logging.info(f"Evidence strings saved into {outputFile}. Exiting.")

if __name__ == '__main__':
    main()