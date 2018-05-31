#!/usr/bin/env python
# -*- coding: utf-8 -*-

from assessment.assessmentTools import Info, read_benchmark
from assessment.GOPrediction import GOPrediction
from assessment.RESULT import result
import helper
import os
import sys
import gc
import pickle as cp
import time

if __name__=='__main__':
    '''
    Main function that takes a predicition and returns calculated values
    '''
    start_time = time.time()
    # Read Config
    obo_path, ic_path, prediction_path, benchmark_directory, results_directory = helper.read_config_MAIN()
    print("Finished reading config ")
    print(time.time() - start_time)
    # Setup workspace
    print('\n Evaluating {}\n'.format(prediction_path))
    # Get predictions
    all_prediction  = GOPrediction()
    prediction_file = open(prediction_path, 'r')
    print("Read in prediction ")
    print(time.time() - start_time)
    # Read in predictions, split by ontology, and save to disk
    all_prediction.read_and_split_and_write(obo_path, prediction_file)
    # Store values
    author, model, keywords, taxon = all_prediction.author, all_prediction.model, all_prediction.keywords, all_prediction.taxon
    # Clear memory
    del all_prediction
    gc.collect()
    # Print values
    print('AUTHOR:   {}\n'.format(author   ))
    print('MODEL:    {}\n'.format(model    ))
    print('KEYWORDS: {}\n'.format(keywords ))
    print('Species:  {}\n'.format(taxon    ))
    print("Done with prediction / Start making result")
    print(time.time() - start_time)
    # Make Results path
    results_path = results_directory + " " + author
    # Make directory 
    helper.mkdir_results(results_path)
    # Make a result object for storing output
    r = result(results_path)
    # Populate result
    r.info(author, model, keywords, taxon, prediction_path)
    print("Done making result")
    print(time.time() - start_time)
    # For each ontology
    for ontology in ['bpo','cco','mfo']:
        print("Starting {}".format(ontology))
        print(time.time() - start_time)
        path = os.path.splitext(prediction_file.name)[0] + '_' + ontology.upper() + '.txt'
        # Grab calculated IC for ontology
        ic_map = cp.load(open("{}ia_{}.map".format(ic_path, ontology.upper()), "rb"))
        print('Ontology: {}\n'.format(ontology))
        # For each type (NK, LK)
        for Type in ['type1','type2']:
            print("Starting {}, reading benchmark".format(Type))
            print(time.time() - start_time)
            print('benchmark type:{}\n' .format(helper.typeConverter(Type)))
            benchmark, obocountDict = read_benchmark(ontology, helper.taxon_name_converter(taxon), Type, benchmark_directory, obo_path)
            if benchmark == None:
                sys.stderr.write('No benchmark is available for the input species and type')
            # Create information object for passing necessary information to subroutines    
            i = Info(benchmark, path, obocountDict[ontology], ic_map, results_path)
            i.writeBenchwithIC(ontology, Type)
            # Check for success
            if i.exist:
                # For each mode
                for mode in ['partial', 'full']:
                    print ('Mode: {}\n'.format(mode))
                    r.coverage = i.coverage()
                    
                    print("Starting FMAX")
                    print(time.time() - start_time)
                    r_temp = r
                    #print ('FMAX')
                    fm = i.check("FMAX", ontology, Type, mode)
                    r_temp.update(fm[0], fm[1], fm[2], fm[3], fm[4])
                    r_temp.writeOut(ontology, Type, mode, "FMAX")
                    
                    print("Starting WFMAX")
                    print(time.time() - start_time)
                    r_temp = r
                    #print ('WFMAX')
                    wfm = i.check("WFMAX", ontology, Type, mode)
                    r_temp.update(wfm[0], wfm[1], wfm[2], wfm[3], wfm[4])
                    r_temp.writeOut(ontology, Type, mode, "WFMAX")
                    
                    print("Starting SMIN")
                    print(time.time() - start_time)
                    r_temp = r
                    #print ('SMIN')                    
                    sm = i.check("SMIN", ontology, Type, mode)
                    r_temp.update(sm[0], sm[1], sm[2], sm[3], sm[4])
                    r_temp.writeOut(ontology, Type, mode, "SMIN")
                    
                    print("Starting NSMIN")
                    print(time.time() - start_time)
                    r_temp = r
                    #print ('NSMIN')
                    nsm = i.check("NSMIN", ontology, Type, mode)
                    r_temp.update(nsm[0], nsm[1], nsm[2], nsm[3], nsm[4])
                    r_temp.writeOut(ontology, Type, mode, "NSMIN")
                    
                    print("Starting AUC")
                    print(time.time() - start_time)
                    r_temp = r
                    print ('AUC')
                    #auc = i.check("AUC", ontology, Type, mode)
                    #r.update("AUC", auc[0], auc[1], auc[2], auc[3], auc[4])
                    #r.writeOut(ontology, Type, mode, "AUC")
               
            del i
            gc.collect()
