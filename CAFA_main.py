#!/usr/bin/env python
# -*- coding: utf-8 -*-


from assessment.GOPrediction import GOPrediction
from assessment.tools import Data, assessMetrics, readOBO, vprint
import helper
import os
import gc

if __name__=='__main__':
    '''
    Main function that takes a predicition and returns calculated values
    '''
    # Read Config
    obo_path, ic_path, prediction_path, benchmark_directory, results_directory = helper.read_config_MAIN()
    # Setup workspace
    vprint('\n Evaluating {}\n'.format(prediction_path), 1 )
    
    ###################################### Prediction IN #################################3
    # Get predictions
    all_prediction  = GOPrediction()
    prediction_file = open(prediction_path, 'r')
    # Read in predictions, split by ontology, and save to disk
    all_prediction.read_and_split_and_write(obo_path, prediction_file)
    # Store values
    author, model, keywords, taxon = all_prediction.author, all_prediction.model, all_prediction.keywords, all_prediction.taxon
    # Make Results path
    results_path = results_directory + " " + author    
    # Memory Management    
    del all_prediction
    gc.collect()
    
    #Read OBO file
    obocounts = readOBO(obo_path)
    vprint(obocounts, 9)
    # RUN ON ALL GROUPS
    # For each ontology
    for ontology in ['BPO','CCO','MFO']:
        path = os.path.splitext(prediction_file.name)[0] + '_' + ontology + '.txt'
        vprint(path, 9)
        # For each type (NK, LK)
        for Type in ['type1','type2']:
            data = Data(ontology, Type)
            data.setTeamInfo(author, model, keywords, taxon)
            data.setPaths(results_path, prediction_path, obo_path, path)
            data.setIC(ic_path)
            data.setBenchmark(taxon, obocounts, benchmark_directory)
            data.setPrediction()
            data.buildLists()
            vprint("Done building benchmark / prep work)", 9)
            # For each mode
            for mode in ['partial', 'full']:
                vprint("Mode: {}".format(mode), 9)
                # RUN METRICS
                r = assessMetrics(data, mode)
                #r.toFile()    
                vprint("Done", 8)
