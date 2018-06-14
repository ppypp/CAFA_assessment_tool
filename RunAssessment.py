#!/usr/bin/env python
# -*- coding: utf-8 -*-

from assessment.Fmetric import FMAX
from assessment.WeightedFmetric import WFMAX
from assessment.Smetric import SMIN
from assessment.NormalizedSmetric import NSMIN
from assessment.GOPrediction import GOPrediction
from assessment.Tools import Info, readOBO, vprint
import helper
import os
import gc
import pickle as cp

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
    for ontology in ['CCO']:
    #for ontology in ['BPO','CCO','MFO']:
        vprint(ontology,1)
        
        path = os.path.splitext(prediction_file.name)[0] + '_' + ontology + '.txt'
        vprint(path, 9)
        
        # For each type (NK, LK)
        for Type in ['type1']:
        #for Type in ['type1','type2']:
            vprint(Type,1)
            info = Info(ontology, Type)
            info.setTeamInfo(author, model, keywords, taxon)
            info.setPaths(results_path, prediction_path, obo_path, path)
            
            cp.dump(info, open("TeamInfo.info","wb"))
            vprint("Saved TeamInfo Info",1)
            #info = cp.load(open("TeamInfo.info", "rb"))
            
            info.setIC(ic_path)
            
            cp.dump(info, open("IC.info","wb"))
            vprint("Saved IC Info",1)
            #info = cp.load( open("IC.info","rb"))
            
            info.setBenchmark(obocounts, benchmark_directory)
            
            cp.dump(info, open("Benchmark.info","wb"))
            vprint("Saved Benchmark Info",1)
            #info = cp.load( open("Benchmark.info","rb"))
            
            info.setPrediction()
            
            cp.dump(info, open("Prediction.info","wb"))
            vprint("Saved Prediction Info",1)
            #info = cp.load(open("Prediction.info","rb"))
            
            info.setInfo()
            
            cp.dump(info, open("Info.info","wb"))
            vprint("Saved Final Info",1)
            #info = cp.load(open("Info.info","rb"))
            
            vprint("Done  prep work)", 9)
            # For each mode
            for mode in ['partial', 'full']:
                info.setMode(mode)
                vprint("Mode: {}".format(mode), 9)
                # RUN METRICS
                #print("FMAX")
                #print(FMAX(info))
                #print("WFMAX")
                #print(WFMAX(info))
                #print("SMIN")
                #print(SMIN(info))
                #print("NSMIN")
                #print(NSMIN(info))
                vprint("Done", 8)
