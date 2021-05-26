import argparse
import logging
import sys
import time

import xml.etree.ElementTree as ET
from pyspark.conf import SparkConf
from pyspark.sql import SparkSession, Row
from pyspark.sql.types import StringType
import pyspark.sql.functions as F

from ontoma import OnToma

# The rest of the types are assigned to -> germline for allele origins
EXCLUDED_ASSOCIATIONTYPES = [
    "Major susceptibility factor in",
    "Part of a fusion gene in",
    "Candidate gene tested in",
    "Role in the phenotype of",
    "Biomarker tested in",
    "Disease-causing somatic mutation(s) in"
]

# Assigning variantFunctionalConsequenceId:
CONSEQUENCE_MAP = {
    "Disease-causing germline mutation(s) (loss of function) in": "SO_0002054",
    "Disease-causing germline mutation(s) in": None,
    "Modifying germline mutation in": None,
    "Disease-causing germline mutation(s) (gain of function) in": "SO_0002053"
}

class ontoma_efo_lookup():
    '''
    Simple class to map orphanet ids to efo
    '''
    def __init__(self):
        self.otmap = OnToma()

    def get_mapping(self, disease_lable, disease_id):

        mappings = self.query_ontoma(disease_id)

        if mappings and 'EFO' in mappings['source']:
            return mappings['term'].split('/')[-1]
        else:
            mappings = self.query_ontoma(disease_id)

        if mappings and 'EFO' in mappings['source']:
            return mappings['term'].split('/')[-1]
        else:
            return None

    def query_ontoma(self, term):

        try:
            mappings = self.otmap.find_term(term, verbose=True)
        except Exception:
            time.sleep(3)
            mappings = self.otmap.find_term(term, verbose=True)

        return mappings

def check_data(orphanet_df) -> None:
    '''
    Function to generate basic stats about the data to the log

    Args:
        orphanet_df (pd.DataFrame): parsed orphanet data with mapped diseases
    '''

    logging.info(f"Number of evidence: {len(orphanet_df)}")
    logging.info(f"Number of unique diseases: {len(orphanet_df.diseaseFromSourceId.unique)}")
    logging.info(f"Number of unique targets: {len(orphanet_df.targetFromSource.unique())}")

    genes_notmapped = orphanet_df.loc[lambda df: df['targetFromSourceIds'].apply(lambda x: len(x) == 0)]
    logging.info(f"Number of targets without ensembl gene ID: {len(genes_notmapped)}")

    disease_unmapped = orphanet_df.loc[df['diseaseFromSourceMappedId'].notnull()]
    logging.info(f"Number of diseases without EFO mapping: {len(disease_unmapped)}")

    # Number of associations


def parserOrphanetXml(orphanet_file: str) -> list:
    '''
    Function to parse Orphanet xml dump and return the parsed
    data as a pandas dataframe.

    Args:
        orphanet_file (str): Orphanet XML filename

    Returns:
        parsed data as a list of dictionary
    '''

    # Reading + validating xml:
    tree = ET.parse(orphanet_file)
    assert isinstance(tree, ET.ElementTree)

    root = tree.getroot()
    assert isinstance(root, ET.Element)

    # Checking if the basic nodes are in the xml structure:

    logging.info(f"There are {root.find('DisorderList').get('count')} disease in the Orphanet xml file.")

    orphanet_disorders = []

    for disorder in root.find('DisorderList').findall('Disorder'):

        # Extracting disease information:
        parsed_disorder = {
            "diseaseFromSource": disorder.find('Name').text,
            "diseaseFromSourceId": 'Orphanet_' + disorder.find('OrphaCode').text,
            "type": disorder.find('DisorderType/Name').text,
        }

        # One disease might be mapped to multiple genes:
        for association in disorder.find('DisorderGeneAssociationList'):

            # For each mapped genes, an evidence is created:
            evidence = parsed_disorder.copy()

            # Not all gene/disease association is backed up by publication:
            try:
                evidence['literature'] = [
                    pmid.replace('[PMID]', '') for pmid in association.find('SourceOfValidation').text.split('_') if '[PMID]' in pmid
                ]
            except AttributeError:
                evidence['literature'] = None

            evidence['associationType'] = association.find('DisorderGeneAssociationType/Name').text
            evidence['confidence'] = association.find('DisorderGeneAssociationStatus/Name').text

            # Parse gene name and id - going for Ensembl gene id only:
            gene = association.find('Gene')
            evidence['targetFromSource'] = gene.find('Name').text

            # Extracting ensembl gene id from cross references:
            ensembl_gene_id = [
                xref.find('Reference').text for xref in gene.find('ExternalReferenceList') if 'ENSG' in xref.find('Reference').text
            ]
            evidence['targetFromSourceId'] = ensembl_gene_id[0] if len(ensembl_gene_id) > 0 else None

            # Collect evidence:
            orphanet_disorders.append(evidence)

    return orphanet_disorders


def main(input_file: str, output_file: str, local: bool = False) -> None:

    # Initialize spark session
    if local:
        sparkConf = (
            SparkConf()
            .set('spark.driver.memory', '15g')
            .set('spark.executor.memory', '15g')
            .set('spark.driver.maxResultSize', '0')
            .set('spark.debug.maxToStringFields', '2000')
            .set('spark.sql.execution.arrow.maxRecordsPerBatch', '500000')
        )
        spark = (
            SparkSession.builder
            .config(conf=sparkConf)
            .master('local[*]')
            .getOrCreate()
        )
    else:
        sparkConf = (
            SparkConf()
            .set('spark.driver.maxResultSize', '0')
            .set('spark.debug.maxToStringFields', '2000')
            .set('spark.sql.execution.arrow.maxRecordsPerBatch', '500000')
        )
        spark = (
            SparkSession.builder
            .config(conf=sparkConf)
            .getOrCreate()
        )

    # Initialize mapping object:
    ol_obj = ontoma_efo_lookup()
    ont_udf = F.udf(ol_obj.get_mapping, StringType())

    # UDF to map association type to sequence ontology:
    so_map = F.udf(lambda x: CONSEQUENCE_MAP[x], StringType())

    # Parsing xml file:s
    orphanet_disorders = parserOrphanetXml(input_file)

    # Crete a spark dataframe from the parsed data:
    orphanet_df = (
        spark.createDataFrame(Row(**x) for x in orphanet_disorders)
        .filter(
            ~F.col('associationType').isin(EXCLUDED_ASSOCIATIONTYPES)
        )
        .withColumn('dataSourceId', F.lit('orphanet'))
        .withColumn('datatypeId', F.lit('genetic_association'))
        .withColumn('AlleleOrigins', F.split(F.lit('germline'), "_"))
        .withColumn('variantFunctionalConsequenceId', so_map(F.col('associationType')))
        .drop(*['associationType', 'type'])
    )

    # Map diseases:
    orphanet_mapping = (
        orphanet_df
        .select('diseaseFromSourceId', 'diseaseFromSource')
        .distinct()
        .withColumn('diseaseFromSourceMappedId', ont_udf(F.col('diseaseFromSourceId'), F.col('diseaseFromSource')))
        .drop('diseaseFromSource')
    )

    # Join mapping:
    orphanet_df = (
        orphanet_df
        .join(orphanet_mapping, on='diseaseFromSourceId', how='right')
    )

    # Save data:
    (
        orphanet_df
        .write.format('json').mode('overwrite').option('compression', 'gzip')
        .save(output_file)
    )


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Parse Orphanet gene-disease annotation downloaded from http://www.orphadata.org/data/xml/en_product6.xml')
    parser.add_argument('--input_file', help='Xml file containing target/disease associations.', type=str)
    parser.add_argument('--output_file', help='Name of the gzipped, JSON evidence file.', type=str)
    parser.add_argument('--logFile', help='Destination of the logs generated by this script.', type=str, required=False)
    parser.add_argument(
        '--local', action='store_true', required=False, default=False,
        help='Flag to indicate if the script is executed locally or on the cluster'
    )

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    log_file = args.logFile
    is_local = args.local

    # Initialize logging:
    logging.basicConfig(
        level=logging.INFO,
        format='%(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if log_file:
        logging.config.fileConfig(filename=log_file)
    else:
        logging.StreamHandler(sys.stderr)

    main(input_file, output_file, is_local)
