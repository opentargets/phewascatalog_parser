#!/usr/bin/env python

import argparse
import sys
from pyspark import *
from pyspark.sql import *
from pyspark.sql.types import *
from pyspark.sql.functions import *
from pyspark.sql.window import Window
import logging


# example run local, make sure you dont pass many parquet files as input to not run out of mem
# also, double check you put the ram size properly in the local settings in the SparkConf map
# python modules/EPMC.py --local \
#        --cooccurrenceFile /home/mkarmona/src/opentargets/data/platform/epmc/epmc-cooccurrences/part-0000\* \
#        --outputFile test_out


def main():

    ##
    ## Parsing parameters:
    ##
    parser = argparse.ArgumentParser(description='This script generates target/disease evidence strings from ePMC cooccurrence files.')
    parser.add_argument('--cooccurrenceFile', help='Partioned parquet file with the ePMC cooccurrences', type=str, required=True)
    parser.add_argument('--outputFile', help='Resulting evidence file saved as compressed JSON.', type=str, required=True)
    parser.add_argument('--logFile', help='Destination of the logs generated by this script.', type=str, required=False)
    parser.add_argument('--local', help='Destination of the logs generated by this script.', action='store_true', required=False, default=False)
    args = parser.parse_args()

    # extract parameters:  
    cooccurrenceFile = args.cooccurrenceFile


    # Initialize logger based on the provided logfile. 
    # If no logfile is specified, logs are written to stderr 
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    if args.logFile:
        logging.config.fileConfig(filename=args.logFile)
    else:
        logging.StreamHandler(sys.stderr)

    # Parse output file:
    out_file = args.outputFile

    ##
    ## Initialize spark session
    ##
    global spark
    local = True

    if args.local:
        sparkConf = (SparkConf()
                 .set("spark.driver.memory", "10g")
                 .set("spark.executor.memory", "10g")
                 .set("spark.driver.maxResultSize", "0")
                 .set("spark.debug.maxToStringFields", "2000")
                 .set("spark.sql.execution.arrow.maxRecordsPerBatch", "500000")
                 )
        spark = (
            SparkSession.builder
                .config(conf=sparkConf)
                .master('local[*]')
                .getOrCreate()
        )
    else:
        sparkConf = (SparkConf()
                 .set("spark.driver.maxResultSize", "0")
                 .set("spark.debug.maxToStringFields", "2000")
                 .set("spark.sql.execution.arrow.maxRecordsPerBatch", "500000")
                 )
        spark = (
            SparkSession.builder
                .config(conf=sparkConf)
                .getOrCreate()
        )
    logging.info(f'Spark version: {spark.version}')

    ##
    ## Log parameters:
    ##
    logging.info(f'Cooccurrence file: {cooccurrenceFile}')
    logging.info(f'Output file: {out_file}')
    logging.info(f'Generating evidence:')

    ##
    ## Load/filter datasets
    ##
    partitionKeys = ['pmid', 'targetFromSourceId', 'diseaseFromSourceMappedId']
    w = Window.partitionBy(*partitionKeys)
    (
        # Reading file:
        spark.read.parquet(cooccurrenceFile)

        # Filtering for diases/target cooccurrences:
        .filter((col('type') == "GP-DS") & (col('isMapped') == True))

        # Renaming columns:
        .withColumnRenamed("keywordId1", "targetFromSourceId")
        .withColumnRenamed("keywordId2", "diseaseFromSourceMappedId")
        .withColumnRenamed("label1", "targetFromSource")
        .withColumnRenamed("label2", "diseaseFromSource")

            # collect sets of field values per window aggregation in w with keys partitionKeys
        .withColumn('textMiningSentences', collect_list(
                struct(
                    col("text"),
                    col('start1').alias('tStart'),
                    col("end1").alias('tEnd'),
                    col('start2').alias('dStart'),
                    col("end2").alias('dEnd'), 
                    col('section')
                )).over(w)
            )
        .withColumn("literature", collect_set(col('pmid')).over(w))
        .withColumnRenamed("label1", "targetFromSource")
        .withColumnRenamed("label2", "diseaseFromSource")
        .dropDuplicates(partitionKeys)

        # Adding linteral columns:
        .withColumn('datasourceId',lit('europepmc'))
        .withColumn('datatypeId',lit('literature'))

        # Reorder columns:
        .select(["datasourceId", "datatypeId", "targetFromSource", "targetFromSourceId",
                "diseaseFromSource","diseaseFromSourceMappedId","literature","textMiningSentences"])

        # Save output:
        .write.format('json').mode('overwrite').option('compression', 'gzip').save(args.outputFile)
    )

    return 0


if __name__ == '__main__':

    main()

