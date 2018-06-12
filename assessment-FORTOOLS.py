#!/usr/bin/env python
# -*- coding: utf-8 -*-

from assessment.assessmentTools import Info, read_benchmark
from assessment.GOPrediction import GOPrediction
from assessment.tools import result
import helper
import os
import sys
import gc
import pickle as cp
import tools
from assessment.tools import Data

if __name__=='__main__':
    '''
    Main function that takes a predicition and returns calculated values
    '''
    # Read Config
    obo_path, ic_path, prediction_path, benchmark_directory, results_directory = helper.read_config_MAIN()

        
    # Setup workspace
    print('\n Evaluating {}\n'.format(prediction_path))
    # Get predictions
    all_prediction  = GOPrediction()
    prediction_file = open(prediction_path, 'r')
    # Read in predictions, split by ontology, and save to disk
    all_prediction.read_and_split_and_write(obo_path, prediction_file)
    # Store values
    author, model, keywords, taxon = all_prediction.author, all_prediction.model, all_prediction.keywords, all_prediction.taxon

    # Make Results path
    results_path = results_directory + " " + author
    # Make a result object for storing output
    # Populate result
    r.teamInfo(author, model, keywords, taxon, prediction_path)
    
    
    
    
    # For each ontology
    for ontology in ['bpo','cco','mfo']:
        path = os.path.splitext(prediction_file.name)[0] + '_' + ontology.upper() + '.txt'
        # For each type (NK, LK)
        for Type in ['type1','type2']:
            data = Data(ontology, Type)
            data.setPaths(results_path, prediction_path, obo_path, path)
            data.setIC(ic_path)
            data.setBenchmark(taxon, benchmark_directory)
            data.setPrediction()
            # For each mode
            for mode in ['partial', 'full']:
                # RUN METRICS
                r = tools.assessMetrics(data)
                
            del i
            gc.collect()
