import ontoma
import logging
import argparse
import json
import gzip

import pandas as pd

class ClinGen():
    def __init__(self):

        self.evidence_strings = list()
        self.unmapped_diseases = set()

        # Create formatter
        formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

        # Create OnToma object
        self.ontoma = ontoma.interface.OnToma()

    def process_gene_validity_curations(self, in_filename, out_filename):

        self.generate_evidence_strings(in_filename)

        # Save results to file
        self.write_evidence_strings(out_filename)

    def generate_evidence_strings(self, filename):

        # When reading csv file skip header lines that don't contain column names
        gene_validity_curation_df = pd.read_csv(filename, skiprows=[0, 1, 2, 3, 5], quotechar='"')

        for index, row in gene_validity_curation_df.iterrows():
            logging.info('{} - {}'.format(row['GENE SYMBOL'], row['DISEASE LABEL']))

            disease_name = row['DISEASE LABEL']
            disease_id = row['DISEASE ID (MONDO)']
            efo_mappings = None

            evidence = {
                'datasourceId': 'clingen',
                'datatypeId': 'genetic_literature',
                'targetFromSourceId': row['GENE SYMBOL'].rstrip(),
                'diseaseFromSource': disease_name,
                'diseaseFromSourceId': disease_id,
                'allelicRequirements': [row['MOI']],
                'confidence': row['CLASSIFICATION'],
                'studyId': row['GCEP'],
                'urls': [{'url': row['ONLINE REPORT']}]
            }

            # Looking up disease in Ontoma using disease id and label:
            if self.ontoma.get_efo_label(disease_id):
                disease_label = self.ontoma.get_efo_label(disease_id)
                # Create list of single disease to mimic what is returned by next step
                efo_mappings = [{'id': disease_id, 'name': disease_label}]

            elif self.ontoma.get_efo_from_xref(disease_id):
                efo_mappings = self.ontoma.get_efo_from_xref(disease_id)
                logging.info('{} mapped to {} EFO ids based on xrefs.'.format(disease_id, len(efo_mappings)))

            else:
                # Search disease label using OnToma and accept perfect matches
                ontoma_mapping = self.ontoma.find_term(disease_name, verbose=True)
                if ontoma_mapping:
                    if ontoma_mapping['action'] is None:
                        efo_mappings = [{'id': ontoma_mapping['term'], 'name': ontoma_mapping['label']}]
                    else:
                        # OnToma fuzzy match ignored
                        logging.info('Fuzzy match from OnToma ignored. Request EFO team to import {} - {}'.format(disease_name, disease_id))
                        # Record the unmapped disease
                        self.unmapped_diseases.add((disease_id, disease_name))
                else:
                    # MONDO id could not be found in EFO. Log it and continue
                    logging.info('{} - {} could not be mapped to any EFO id. Skipping it, it should be checked with the EFO team'.format(disease_name, disease_id))
                    # Record the unmapped disease
                    self.unmapped_diseases.add((disease_id, disease_name))

            # Generating evidence for all mapped efo:
            if efo_mappings:
                for efo_mapping in efo_mappings:
                    evidence_with_efo = evidence.copy()
                    evidence_with_efo['diseaseFromSourceMappedId'] = ontoma.interface.make_uri(efo_mapping['id']).split('/')[-1]
                    self.evidence_strings.append(evidence_with_efo)
            else:
                self.evidence_strings.append(evidence)

        if len(self.unmapped_diseases) > 0:
            logging.info(f'There are {len(self.unmapped_diseases)} unmapped diseases.')
            unmapped_diseases_string = "\n- ".join([x[1] for x in self.unmapped_diseases])
            logging.info(f'Unmapped diseases: \n- {unmapped_diseases_string}')

    def write_evidence_strings(self, filename):
        logging.info('Writing ClinGen evidence strings to %s', filename)
        with gzip.open(filename, 'wt') as tp_file:
            for evidence_string in self.evidence_strings:
                json.dump(evidence_string, tp_file)
                tp_file.write('\n')


def main(infile, outfile):

    clingen = ClinGen()
    clingen.process_gene_validity_curations(infile, outfile)


if __name__ == "__main__":

    # Parse CLI arguments
    parser = argparse.ArgumentParser(description='Parse ClinGen gene-disease associations from Gene Validity Curations')
    parser.add_argument('-i', '--input_file',
                        help='Name of csv file downloaded from https://search.clinicalgenome.org/kb/gene-validity',
                        type=str, required=True)
    parser.add_argument('-o', '--output_file',
                        help='Name of evidence file',
                        type=str, required=True)
    parser.add_argument('-l', '--log_file', type=str,
                        help='Optional filename to redirect the logs into.')

    args = parser.parse_args()

    # Get parameters
    infile = args.input_file
    outfile = args.output_file
    log_file = args.log_file

    # Initialize the logger based on the provided log file. If no log file is specified, logs are written to STDERR.
    logging_config = {
        'level': logging.INFO,
        'format': '%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        'datefmt': '%Y-%m-%d %H:%M:%S',
    }
    if log_file:
        logging_config['filename'] = log_file
    logging.basicConfig(**logging_config)

    main(infile, outfile)
