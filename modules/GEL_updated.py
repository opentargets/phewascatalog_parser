import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from ontoma import OnToma

from settings import Config
from datetime import date
import logging
import requests
import python_jsonschema_objects as pjo


class panelapp_evidence_generator():

    def __init__(self, schema_version=Config.VALIDATED_AGAINST_SCHEMA_VERSION):

        # Logger configuration

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        # create console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        # create formatter
        formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
        # add formatter to ch
        ch.setFormatter(formatter)
        # add ch to logger
        self.logger.addHandler(ch)

        # Get JSON Schema

        self.schema_version = schema_version
        schema_url = "https://raw.githubusercontent.com/opentargets/json_schema/" + self.schema_version + "/draft4_schemas/opentargets.json"
        self.logger.info(f"Loading JSON Schema from {schema_url}")


        # Initialize JSON Schema builder

        try:
            r = requests.get(schema_url)
            r.raise_for_status()
            json_schema = r.json()
            self.builder = pjo.ObjectBuilder(json_schema)
            self.evidence_builder = self.builder.build_classes()
        except requests.exceptions.HTTPError as e:
            self.logger.error('Invalid JSON schema version')
            raise e


        # Create OnToma object

        self.otmap = OnToma()

        # Parse Gene Symbol to EnsEMBL ID

        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()
        self.genes = gene_parser.genes

    def build_publications(self, dataframe):
        '''
        Builds a dataframe with the publications fetched from the PanelApp API
        '''
        panel_lst = dataframe["Panel Id"].value_counts().index.to_list()
        lst = []

        for panel in panel_lst:
            tmp = dataframe[dataframe["Panel Id"] == panel]
            request = publications_from_panel(panel)
            tmp["Publications"] = tmp.apply(lambda X: publication_from_symbol(X.Symbol, request), axis=1)
            lst.append(tmp)
        
        dataframe = pd.concat(lst, ignore_index=True, sort=False)
        return dataframe

    
    def publications_from_panel(self, panel_id):
        '''
        queries the PanelApp API to obtain a list of the publications for every gene within a panel_id
        '''
        try:
            url = f"http://panelapp.genomicsengland.co.uk/api/v1/panels/{panel_id}/"
            res = requests.get(url).json()

            return res["genes"]
        except:
            return None
    
    def publication_from_symbol(self, symbol, response):
        '''

        '''
        for gene in response:
            if gene["gene_data"]["gene_symbol"] == symbol:
                publication = gene["publications"]
                return publication
            continue


    def split_dataframes(self, dataframe):
        '''
        Cleaning of the initial .tsv or .csv file
        '''

        # NaNs in "Phenotypes" column --> Assignment of Pannel Name

        dataframe.loc[dataframe["Phenotypes"].isna(), "Phenotypes"] = dataframe["Panel Name"]

        # Splitting between one and multiple associations 
        # "Phenotypes" column contains the separator ";"

        one_phenotype = dataframe[~dataframe["Phenotypes"].str.contains(";", na=False)].reset_index(drop=True)
        multiple_phenotype = dataframe[dataframe["Phenotypes"].str.contains(";", na=False)].reset_index(drop=True)

        # Among one_phenotype, splitting between those that include an OMIM code and those who don't.

        OMIM_codes = one_phenotype[one_phenotype["Phenotypes"].str.contains('[0-9]{5}', na=False, regex=True)].reset_index(drop=True)
        not_OMIM = one_phenotype[~one_phenotype["Phenotypes"].str.contains('[0-9]{5}', na=False, regex=True)].reset_index(drop=True)

        # Cleaning the OMIM_codes df: 
        #   - cleaned_OMIM: df after the removal of the OMIM codes in the Phenotypes column and the inclusion of these in a Codes column
        #   - Deleting special characters (Code is crashing, solution is a shame - TO FIX)

        cleaned_OMIM = OMIM_codes.copy()
        cleaned_OMIM["Codes"] = cleaned_OMIM["Phenotypes"].str.extract(r"(\d{6})")
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].replace("(\d{6})","",regex=True)

        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("{", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("}", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace(",", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("(", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace(")", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("[", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("]", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("'", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace('"', ''))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace(".", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("-", " "))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("?", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.replace("#", ""))
        cleaned_OMIM["Phenotypes"] = cleaned_OMIM["Phenotypes"].apply(lambda X: X.strip())

        # Exporting dfs to .csv

        #cleaned_OMIM.to_csv("cleaned_OMIM.csv", index=False)
        #not_OMIM.to_csv("not_OMIM.csv", index=False)
        #multiple_phenotype.to_csv("multiple_phenotype.csv", index=False)

        return cleaned_OMIM, not_OMIM, multiple_phenotype

    def ontoma_query(self, iterable, dict_name="ontoma_queries.json"):
        '''
        Queries the OnToma utility to map a phenotype to a disease.
        '''

        dct = dict()

        for e in iterable:
            try:
                tmp = self.otmap.find_term(e, verbose=True)
                if tmp != None:
                    dct[e] = tmp
                else:
                    dct[e] = {
                    'term': None,
                     'label': None,
                     'source': None,
                     'quality': None,
                     'action': None
                    }
            except:
                continue
        
        with open(dict_name, "w") as outfile:  
            json.dump(dct, outfile)

        return dct
    

    def get_evidence_string(self, row):


        # Association level information:
        self.gene_symbol = row["Symbol"]
        self.ensembl_iri = "http://identifiers.org/ensembl/" + self.genes[self.gene_symbol]
        self.mapped_disease = row["OnToma Label"]
        self.mapped_id = row["OnToma Term"]
        self.reported_disease = row["Phenotypes"]
        self.mode_of_inheritance = row["Mode of inheritance"]
        self.evidence_confidence = row["List"]
        self.sources = row["Sources"]
        self.publication = row["Publications"]

        # Panel level information:
        self.panel_id = row["Panel Id"]
        self.panel_name = row["Panel Name"]
        self.panel_version = row["Panel Version"]
        self.panel_type = row["Panel Types"]
        self.super_panel_id = row["Super Panel Id"]
        self.super_panel_name = row["Super Panel Name"]
        self.super_panel_version = row["Super Panel Version"]

        # If the association has no EFO attached:
        if self.mapped_id is None or self.mapped_id == '':
            logging.warning('No EFO id for association row: {}'.format(row.name))
            return None

        target_field = {
                'id' : ensembl_iri,
                'activity' : "http://identifiers.org/cttv.activity/unknown",
                'target_type' : "http://identifiers.org/cttv.target/gene_evidence",
                'target_name' : gene_symbol
        }

        provenance_type = {
                        'database' : {
                            'id' : "Genomics England PanelApp",
                            'version' : date.today().strftime("%d/%m/%Y"),
                            'dbxref' : {
                                'url': "https://panelapp.genomicsengland.co.uk/",
                                'id' : "Genomics England PanelApp",
                                'version' : date.today().strftime("%d/%m/%Y")
                                        }
                                    }
                            }
        resource_score = {
                        'type': "probability",
                        'value': 1
                    }
        
        linkout = [
                {
                'url' : "http://europepmc.org/abstract/MED/" + self.publication,
                'nice_name' : 'Gene Validity Curations: {} - {} report'.format(gene_symbol, disease_name)
                }
            ]

        evidence_field = {
                    'is_associated' : True,
                    'confidence' : 1,
                    'allelic_requirement' : self.mode_of_inheritance,
                    'evidence_codes' : ["http://purl.obolibrary.org/obo/ECO_0000204"],
                    'provenance_type' : provenance_type,
                    'date_asserted' : date.today().strftime("%d/%m/%Y"),
                    'resource_score' : resource_score,
                    'urls' : linkout
                    }

        target = {
                        'id' : self.ensembl_iri,
                        'activity' : "http://identifiers.org/cttv.activity/unknown",
                        'target_type' : "http://identifiers.org/cttv.target/gene_evidence",
                        'target_name' : gene_symbol
                    }

        disease_info = {
                        'id': self.mapped_id,
                        'name' : self.mapped_disease,
                        'source_name': self.mapped_disease
                    }
        
        unique_association_fields = {
                        'panel_name' : self.panel_name,
                        'panel_id' : self.panel_id,
                        'panel_version': self.panel_version,
                        'original_disease': self.reported_disease,
                        'target_id' : self.ensembl_iri,
                        'disease_id' : self.mapped_id,
                    }

        

        try:
            evidence = self.evidence_builder.Opentargets(
                type = "genetic_literature",
                access_level = "public",
                sourceID = "genomics_england",
                evidence = evidence_field,
                target = target_field,
                disease = disease_info,
                unique_association_fields = unique_association_fields,
                validated_against_schema_version = self.schema_version
            )
            self.evidence_strings.append(evidence)
        except:
            self._logger.warning('Evidence generation failed for row: {}'.format(c))
            raise

def main(dataframe):

    # Initialize evidence builder object:
    evidence_builder = panelapp_evidence_generator(json_schema="1.6.9")

    # Filtering with green and ambar evidences

    dataframe = dataframe[((dataframe["List"] == "green") | (dataframe["List"] == "ambar")) & (dataframe["Panel Version"] > 1) & (dataframe["Panel Status"] == "PUBLIC")].reset_index(drop=True)
    dataframe.dropna(axis=1, how="all", inplace=True)


    # Feed the dataframe with publications

    dataframe = evidence_builder.build_publications(dataframe)

    # Splitting and cleaning the dataframe according to the phenotype string

    OMIM_codes, not_OMIM, multiple_phenotype = evidence_builder.split_dataframes(dataframe)

    # Mapping the phenotypes to an EFO code
        # (Working only with OMIM)

    dataframe = OMIM_codes.copy()

    phenotypes = dataframe["Phenotypes"].value_counts().index.to_list()

    mappings_dict = evidence_builder.ontoma_query(phenotypes)
    assert mappings_dict, "Disease mappings completed."

    # Splitting dictionaries between fuzzies and matches

    fuzzy = {}
    match = {}

    for i in mappings_dict.keys():
        try:
            if mappings_dict[i]["quality"] == "fuzzy":
                fuzzy[i] = mappings_dict[i]
            elif mappings_dict[i]["quality"] == "match":
                match[i] = mappings_dict[i]
        except:
            continue

    # Creating the corresponding OnToma result for each phenotype
    # New columns: OnToma Result, OnToma Source, OnToma Term, OnToma Label
    # TO-DO: Refactor this

    for phenotype in match.keys():
        dataframe.loc[dataframe["mesh_heading"] == phenotype, "OnToma Result"] = "match"
        dataframe.loc[dataframe["mesh_heading"] == phenotype, "OnToma Source"] = match[phenotype]["source"]
        dataframe.loc[dataframe["mesh_heading"] == phenotype, "OnToma Term"] = match[phenotype]["term"]
        dataframe.loc[dataframe["mesh_heading"] == phenotype, "OnToma Label"] = match[phenotype]["label"]
        
    for phenotype in fuzzy.keys():
        dataframe.loc[dataframe["mesh_heading"] == phenotype, "OnToma Result"] = "fuzzy"
        dataframe.loc[dataframe["mesh_heading"] == phenotype, "OnToma Source"] = fuzzy[phenotype]["source"]
        dataframe.loc[dataframe["mesh_heading"] == phenotype, "OnToma Term"] = fuzzy[phenotype]["term"]
        dataframe.loc[dataframe["mesh_heading"] == phenotype, "OnToma Label"] = fuzzy[phenotype]["label"]
        

    evidences = dataframe.apply(evidence_builder.get_evidence_string, axis=1)

    return evidences

test = main("genes_toy.csv")

print(test)

