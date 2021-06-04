#!/usr/bin/env python3
"""Evidence parser for the animal model sources from PhenoDigm."""

import argparse
import logging
import os
import pathlib
import shutil
import tempfile
import urllib.request

import pyspark
import pyspark.sql.functions as pf
import requests
from retry import retry


# The tables and their fields to fetch from SOLR. Other tables (not currently used): gene, disease_gene_summary.
IMPC_SOLR_TABLES = {
    # Mouse to human mappings.
    'gene_gene': ('gene_id', 'hgnc_gene_id'),
    'ontology_ontology': ('mp_id', 'hp_id'),
    # Mouse model and disease data.
    'mouse_model': ('model_id', 'model_phenotypes'),
    'disease': ('disease_id', 'disease_phenotypes'),
    'disease_model_summary': ('model_id', 'model_genetic_background', 'model_description', 'disease_id', 'disease_term',
                              'disease_model_avg_norm', 'disease_model_max_norm', 'marker_id'),
    'ontology': ('ontology', 'phenotype_id', 'phenotype_term'),
}


class ImpcSolrRetriever:
    """Retrieve data from the IMPC SOLR API and save the CSV files to the specified location."""

    # Mouse model data from IMPC SOLR.
    IMPC_SOLR_HOST = 'http://www.ebi.ac.uk/mi/impc/solr/phenodigm/select'

    # The largest table is about 7 million records. The one billion limit is used as an arbitrary high number to
    # retrieve all records in one large request, which maximises the performance.
    IMPC_SOLR_BATCH_SIZE = 1000000000
    IMPC_SOLR_TIMEOUT = 3600

    # The decorator ensures that the requests are retried in case of network or server errors.
    @retry(tries=3, delay=5, backoff=1.2, jitter=(1, 3))
    def get_number_of_solr_records(self, data_type):
        params = {'q': '*:*', 'fq': f'type:{data_type}', 'rows': 0}
        response = requests.get(self.IMPC_SOLR_HOST, params=params, timeout=self.IMPC_SOLR_TIMEOUT)
        response.raise_for_status()  # Check for HTTP errors. This will be caught by @retry.
        return response.json()['response']['numFound']

    @retry(tries=3, delay=5, backoff=1.2, jitter=(1, 3))
    def query_solr(self, data_type, start):
        """Request one batch of SOLR records of the specified data type and write it into a temporary file."""
        list_of_columns = [column.split(' > ')[0] for column in IMPC_SOLR_TABLES[data_type]]
        params = {'q': '*:*', 'fq': f'type:{data_type}', 'start': start, 'rows': self.IMPC_SOLR_BATCH_SIZE, 'wt': 'csv',
                  'fl': ','.join(list_of_columns)}
        response = requests.get(self.IMPC_SOLR_HOST, params=params, timeout=self.IMPC_SOLR_TIMEOUT, stream=True)
        response.raise_for_status()
        # Write records as they appear to avoid keeping the entire response in memory.
        with tempfile.NamedTemporaryFile('wt', delete=False) as tmp_file:
            response_lines = response.iter_lines(decode_unicode=True)
            header = next(response_lines)
            if start == 0:  # Only write the header for the first requested batch.
                tmp_file.write(header + '\n')
            number_of_records = 0
            for line in response_lines:
                number_of_records += 1
                tmp_file.write(line + '\n')
            return number_of_records, tmp_file.name

    def fetch_data(self, data_type, output_filename):
        """Fetch all rows of the requested data type to the specified location."""
        total_records = self.get_number_of_solr_records(data_type)
        assert total_records != 0, f'SOLR did not return any data for {data_type}.'
        with open(output_filename, 'wb') as outfile:
            start, total = 0, 0  # Initialise the counters.
            while True:
                number_of_records, tmp_filename = self.query_solr(data_type, start)
                with open(tmp_filename, 'rb') as tmp_file:
                    shutil.copyfileobj(tmp_file, outfile)
                os.remove(tmp_filename)
                # Increment the counters.
                start += self.IMPC_SOLR_BATCH_SIZE
                total += number_of_records
                # Exit when all documents have been retrieved.
                if total == total_records:
                    break


class PhenoDigm:
    """Retrieve the data, load it into Spark, process and write the resulting evidence strings."""

    # Human gene mappings.
    HGNC_DATASET_URI = 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt'
    HGNC_DATASET_FILENAME = 'hgnc_complete_set.txt'

    # Mouse gene mappings.
    MGI_DATASET_URI = 'http://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt'
    MGI_DATASET_FILENAME = 'MGI_Gene_Model_Coord.rpt'

    IMPC_FILENAME = 'impc_solr_{data_type}.csv'

    def __init__(self, logger, cache_dir):
        self.logger = logger
        self.cache_dir = cache_dir
        self.spark = pyspark.sql.SparkSession.builder.appName('phenodigm_parser').getOrCreate()
        self.hgnc_gene_id_to_ensembl_human_gene_id, self.mgi_gene_id_to_ensembl_mouse_gene_id = [None] * 2
        self.mouse_gene_to_human_gene, self.mouse_phenotype_to_human_phenotype = [None] * 2
        self.mouse_model, self.disease, self.disease_model_summary, self.ontology = [None] * 4
        self.evidence = None

    def update_cache(self):
        """Fetch the Ensembl gene ID and SOLR data into the local cache directory."""
        pathlib.Path(self.cache_dir).mkdir(parents=False, exist_ok=True)

        self.logger.info('Fetching human gene ID mappings from HGNC.')
        urllib.request.urlretrieve(self.HGNC_DATASET_URI, os.path.join(self.cache_dir, self.HGNC_DATASET_FILENAME))

        self.logger.info('Fetching mouse gene ID mappings from MGI.')
        urllib.request.urlretrieve(self.MGI_DATASET_URI, os.path.join(self.cache_dir, self.MGI_DATASET_FILENAME))

        self.logger.info('Fetching PhenoDigm data from IMPC SOLR.')
        impc_solr_retriever = ImpcSolrRetriever()
        for data_type in IMPC_SOLR_TABLES:
            self.logger.info(f'Fetching PhenoDigm data type {data_type}.')
            filename = os.path.join(self.cache_dir, self.IMPC_FILENAME.format(data_type=data_type))
            impc_solr_retriever.fetch_data(data_type, filename)

    def load_tsv(self, filename):
        return self.spark.read.csv(os.path.join(self.cache_dir, filename), sep='\t', header=True, nullValue='null')

    def load_solr_csv(self, data_type):
        """Load the CSV from SOLR; rename and select columns as specified."""
        df = self.spark.read.csv(
            os.path.join(self.cache_dir, self.IMPC_FILENAME.format(data_type=data_type)),
            header=True
        )
        column_name_mappings = [column_map.split(' > ') for column_map in IMPC_SOLR_TABLES[data_type]]
        columns_to_rename = {mapping[0]: mapping[1] for mapping in column_name_mappings if len(mapping) == 2}
        new_column_names = [mapping[-1] for mapping in column_name_mappings]
        # Rename columns.
        for old_column_name, new_column_name in columns_to_rename.items():
            df = df.withColumnRenamed(old_column_name, new_column_name)
        # Restrict only to the columns we need.
        return df.select(new_column_names)

    def load_data_from_cache(self):
        """Load the Ensembl gene ID and SOLR data from the downloaded TSV/CSV files into Spark."""
        # Mappings from HGNC/MGI gene IDs to Ensembl gene IDs.
        self.hgnc_gene_id_to_ensembl_human_gene_id = (  # E.g. 'HGNC:5', 'ENSG00000121410'.
            self.load_tsv(self.HGNC_DATASET_FILENAME)
            .withColumnRenamed('hgnc_id', 'hgnc_gene_id')
            .withColumnRenamed('ensembl_gene_id', 'targetFromSourceId')  # Using the final name.
            .select('hgnc_gene_id', 'targetFromSourceId')
        )
        self.mgi_gene_id_to_ensembl_mouse_gene_id = (  # E.g. 'MGI:87853', 'ENSMUSG00000027596'.
            self.load_tsv(self.MGI_DATASET_FILENAME)
            .withColumnRenamed('1. MGI accession id', 'mgi_gene_id')
            .withColumnRenamed('3. marker symbol', 'targetInModel')  # Using the final name.
            .withColumnRenamed('11. Ensembl gene id', 'targetInModelId')  # Using the final name.
            .filter(pf.col('targetInModelId').isNotNull())
            .select('mgi_gene_id', 'targetInModel', 'targetInModelId')
        )

        # Mouse to human gene mappings, e.g. 'MGI:1346074', 'HGNC:4024'.
        self.mouse_gene_to_human_gene = (
            self.load_solr_csv('gene_gene')
            .withColumnRenamed('gene_id', 'mgi_gene_id')
        )
        # Mouse to human phenotype mappings, e.g. 'MP:0000745','HP:0100033'.
        self.mouse_phenotype_to_human_phenotype = self.load_solr_csv('ontology_ontology')

        # Mouse model and disease data.
        # Note that the models are accessioned with the same prefix ('MGI:') as genes, but they are separate entities.
        self.mouse_model = self.load_solr_csv('mouse_model')  # E. g. 'MGI:3800884', ['MP:0001304 cataract'].
        self.disease = self.load_solr_csv('disease')  # E.g. 'OMIM:609258', ['HP:0000545 Myopia'].
        # E. g. 'MGI:2681494', 'C57BL/6JY-smk', 'smk/smk', 'ORPHA:3097', 'Meacham Syndrome', 91.6, 'MGI:98324'.
        self.disease_model_summary = (
            self.load_solr_csv('disease_model_summary')
            .withColumnRenamed('model_genetic_background', 'biologicalModelGeneticBackground')
            .withColumnRenamed('model_description', 'biologicalModelAllelicComposition')
            # In Phenodigm, the scores report the association between diseases and animal models, not genes. The
            # phenotype similarity is computed using an algorithm called OWLSim which expresses the similarity in terms
            # of the Jaccard Index (simJ) or Information Content (IC). Therefore, to compute the score you can take the
            # maximum score of both analyses (disease_model_max_norm) or a combination of them both
            # (disease_model_avg_norm). In the Results and discussion section of the Phenodigm paper, the methods are
            # compared to a number of gold standards. It is concluded that the geometric mean of both analyses is the
            # superior metric and should therefore be used as the score.
            .withColumnRenamed('disease_model_avg_norm', 'resourceScore')
            .withColumnRenamed('marker_id', 'mgi_gene_id')
        )
        self.ontology = (
            self.load_solr_csv('ontology')  # E.g. 'HP', 'HP:0000002', 'Abnormality of body height'.
            .filter((pf.col('ontology') == 'MP') | (pf.col('ontology') == 'HP'))
        )
        assert self.ontology.select('phenotype_id').distinct().count() == self.ontology.count(), \
            f'Encountered multiple names for the same term in the ontology table.'

    def generate_phenodigm_evidence_strings(self, score_cutoff, ontology_mapping):
        """Generate the evidence by renaming, transforming and joining the columns."""
        # Process ontology information to enable MP and HP term lookup based on the ID.
        mp_terms, hp_terms = (
            self.ontology
            .filter(pf.col('ontology') == ontology_name)
            .withColumnRenamed('phenotype_id', f'{ontology_name.lower()}_id')
            .withColumnRenamed('phenotype_term', f'{ontology_name.lower()}_term')
            .select(f'{ontology_name.lower()}_id', f'{ontology_name.lower()}_term')
            for ontology_name in ('MP', 'HP')
        )

        # Split lists of phenotypes in the `mouse_model` and `disease` tables and keep only the ID. For example, one row
        # with 'MP:0001529 abnormal vocalization,MP:0002981 increased liver weight' becomes two rows with 'MP:0001529'
        # and 'MP:0002981'.
        model_mouse_phenotypes = (
            self.mouse_model
            .withColumn('mp_id', pf.expr(r"regexp_extract_all(model_phenotypes, '(MP:\\d+)', 1)"))
            .withColumn('mp_id', pf.explode('mp_id'))
            .select('model_id', 'mp_id')
        )
        disease_human_phenotypes = (
            self.disease
            .withColumn('hp_id', pf.expr(r"regexp_extract_all(disease_phenotypes, '(HP:\\d+)', 1)"))
            .withColumn('hp_id', pf.explode('hp_id'))
            .select('disease_id', 'hp_id')
        )
        # Map mouse model phenotypes into human terms.
        model_human_phenotypes = (
            model_mouse_phenotypes
            .join(self.mouse_phenotype_to_human_phenotype, on='mp_id', how='inner')
            .select('model_id', 'hp_id')
        )

        # We are reporting all mouse phenotypes for a model, regardless of whether they can be mapped into any human
        # disease.
        all_mouse_phenotypes = (
            model_mouse_phenotypes
            .join(mp_terms, on='mp_id', how='inner')
            .groupby('model_id')
            .agg(
                pf.collect_set(pf.struct(
                    pf.col('mp_id').alias('id'),
                    pf.col('mp_term').alias('label')
                )).alias('diseaseModelAssociatedModelPhenotypes')
            )
            .select('model_id', 'diseaseModelAssociatedModelPhenotypes')
        )
        # For human phenotypes, we only want to include the ones which are present in the disease *and* also can be
        # traced back to the model phenotypes through the MP → HP mapping relationship.
        matched_human_phenotypes = (
            # We start with all possible pairs of model-disease associations.
            self.disease_model_summary.select('model_id', 'disease_id')
            # Add all disease phenotypes. Now we have: model_id, disease_id, hp_id (from disease).
            .join(disease_human_phenotypes, on='disease_id', how='inner')
            # Only keep the phenotypes which also appear in the mouse model (after mapping).
            .join(model_human_phenotypes, on=['model_id', 'hp_id'], how='inner')
            # Add ontology terms in addition to IDs. Now we have: model_id, disease_id, hp_id, hp_term.
            .join(hp_terms, on='hp_id', how='inner')
            .groupby('model_id', 'disease_id')
            .agg(
                pf.collect_set(pf.struct(
                    pf.col('hp_id').alias('id'),
                    pf.col('hp_term').alias('label')
                )).alias('diseaseModelAssociatedHumanPhenotypes')
            )
            .select('model_id', 'disease_id', 'diseaseModelAssociatedHumanPhenotypes')
        )

        self.evidence = (
            # This table contains all unique (model_id, disease_id) associations which form the base of the evidence
            # strings.
            self.disease_model_summary

            # Filter out the associations with a low score.
            .filter(~(pf.col('resourceScore') < score_cutoff))

            # Add the mouse gene mapping information. The mappings are not necessarily one to one, because a single MGI
            # can map to multiple Ensembl mouse genes. When this happens, join will handle the necessary explosions, and
            # a single row from the original table will generate multiple evidence strings. This adds the fields
            # `targetInModel` and `targetInModelId`.
            .join(self.mgi_gene_id_to_ensembl_mouse_gene_id, on='mgi_gene_id', how='inner')
            # Add the human gene mapping information. This is added in two stages: MGI → HGNC → Ensembl human gene.
            # Similarly to mouse gene mappings, at each stage there is a possibility of a row explosion.
            .join(self.mouse_gene_to_human_gene, on='mgi_gene_id', how='inner')
            .join(self.hgnc_gene_id_to_ensembl_human_gene_id, on='hgnc_gene_id', how='inner')  # `targetFromSourceId`.
            .drop('mgi_gene_id', 'hgnc_gene_id')

            # Add all mouse phenotypes of the model → `diseaseModelAssociatedModelPhenotypes`.
            .join(all_mouse_phenotypes, on='model_id', how='left')
            # Add the matched model/disease human phenotypes → 'diseaseModelAssociatedHumanPhenotypes`.
            .join(matched_human_phenotypes, on=['model_id', 'disease_id'], how='left')

            # Model ID adjustments. First, strip the trailing modifiers, where present. The original ID, used for table
            # joins, may look like 'MGI:6274930#hom#early', where the first part is the allele ID and the second
            # specifies the zygotic state. There can be several models for the same allele ID with different phenotypes.
            # However, this information is also duplicated in `biologicalModelGeneticBackground` (for example:
            # 'C57BL/6NCrl,Ubl7<em1(IMPC)Tcp> hom early'), so in this field we strip those modifiers.
            .withColumn(
                'biologicalModelId',
                pf.split(pf.col('model_id'), '#').getItem(0)
            )
            .drop('model_id')
            # Second, we only want to output the model names from the MGI namespace. An example of something we *don't*
            # want is 'NOT-RELEASED-025eb4a791'. This will be converted to null.
            .withColumn(
                'biologicalModelId',
                pf.when(pf.col('biologicalModelId').rlike(r'^MGI:\d+$'), pf.col('biologicalModelId'))
            )

            # Rename the disease data columns.
            .withColumnRenamed('disease_id', 'diseaseFromSourceId')
            .withColumnRenamed('disease_term', 'diseaseFromSource')
            # Add constant value columns.
            .withColumn('datasourceId', pf.lit('phenodigm'))
            .withColumn('datatypeId', pf.lit('animal_model'))

            # Ensure stable column order.
            .select('biologicalModelAllelicComposition', 'biologicalModelGeneticBackground', 'biologicalModelId',
                    'datasourceId', 'datatypeId', 'diseaseFromSource', 'diseaseFromSourceId',
                    'diseaseModelAssociatedHumanPhenotypes', 'diseaseModelAssociatedModelPhenotypes', 'resourceScore',
                    'targetFromSourceId', 'targetInModel', 'targetInModelId')
        )

        # Merge ontology mapping information
        ontology_mapping = self.spark.read.csv(ontology_mapping, sep='\t', header=True, nullValue='null')
        self.evidence = (
            self.evidence
            .join(ontology_mapping, on='diseaseFromSourceId', how='left')
        )

    def write_evidence_strings(self, evidence_strings_filename):
        """Dump the Spark evidence dataframe as a compressed JSON file. The order of the evidence strings is not
        maintained, and they are returned in random order as collected by Spark."""
        with tempfile.TemporaryDirectory() as tmp_dir_name:
            (
                self.evidence.coalesce(1).write.format('json').mode('overwrite')
                .option('compression', 'org.apache.hadoop.io.compress.GzipCodec').save(tmp_dir_name)
            )
            json_chunks = [f for f in os.listdir(tmp_dir_name) if f.endswith('.json.gz')]
            assert len(json_chunks) == 1, f'Expected one JSON file, but found {len(json_chunks)}.'
            os.rename(os.path.join(tmp_dir_name, json_chunks[0]), evidence_strings_filename)


def main(cache_dir, output, score_cutoff, use_cached=False, log_file=None, ontology_mapping=None):
    # Initialize the logger based on the provided log file. If no log file is specified, logs are written to STDERR.
    logging_config = {
        'level': logging.INFO,
        'format': '%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        'datefmt': '%Y-%m-%d %H:%M:%S',
    }
    if log_file:
        logging_config['filename'] = log_file
    logging.basicConfig(**logging_config)

    # Process the data.
    phenodigm = PhenoDigm(logging, cache_dir)
    if not use_cached:
        logging.info('Update the HGNC/MGI/SOLR cache.')
        phenodigm.update_cache()

    logging.info('Load gene mappings and SOLR data from local cache.')
    phenodigm.load_data_from_cache()

    logging.info('Build the evidence strings.')
    phenodigm.generate_phenodigm_evidence_strings(score_cutoff, ontology_mapping)

    logging.info('Collect and write the evidence strings.')
    phenodigm.write_evidence_strings(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    req = parser.add_argument_group('required arguments')
    req.add_argument('--cache-dir', help='Directory to store the HGNC/MGI/SOLR cache files in.', required=True)
    req.add_argument('--output', help='Name of the json.gz file to output the evidence strings into.', required=True)
    parser.add_argument('--score-cutoff', help=(
        'Discard model-disease associations with the `disease_model_avg_norm` score less than this value. The score '
        'ranges from 0 to 100.'
    ), type=float, default=0.0)
    parser.add_argument('--use-cached', help='Use the existing cache and do not update it.', action='store_true')
    parser.add_argument('--log-file', help='Optional filename to redirect the logs into.')

    parser.add_argument('--ontology-mapping')

    args = parser.parse_args()
    main(args.cache_dir, args.output, args.score_cutoff, args.use_cached, args.log_file, args.ontology_mapping)
