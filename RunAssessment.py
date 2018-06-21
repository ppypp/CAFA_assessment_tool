#!/usr/bin/env python
# -*- coding: utf-8 -*-

from assessment_new.Fmetric import FMAX
from assessment_new.WeightedFmetric import WFMAX
from assessment_new.Smetric import SMIN
from assessment_new.NormalizedSmetric import NSMIN
from assessment_new.GOPrediction import GOPrediction
from assessment_new.Tools import Info, readOBO, vprint, vwrite, clear, getTime
import helper
import os
import gc
import pickle as cp

if __name__=='__main__':
    '''
    Main function that takes a predicition and returns calculated values
    '''
    start_time = getTime(0)
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
    
    
     ##################################################
    # Speed up if rerunning
    cp.dump(all_prediction, open("Temp/Prediction.all","wb"))
    #all_prediction = cp.load( open("Temp/Prediction.all","rb"))
    ##################################################
    
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
    #for ontology in ['MFO']:
    for ontology in ['BPO','CCO','MFO']:
        vprint(ontology, 1)
        getTime(start_time)
        path = os.path.splitext(prediction_file.name)[0] + '_' + ontology + '.txt'
        vprint(path, 9)
        
        # For each type (NK, LK)
        #for Type in ['type1']:
        for Type in ['type1','type2']:
            vprint(Type, 1)
            getTime(start_time)
            info = Info(ontology, Type)
            info.setTeamInfo(author, model, keywords, taxon)
            info.setPaths(results_path, prediction_path, obo_path, path)
            info.setTime(start_time)
            
            cp.dump(info, open("Temp/TeamInfo.info","wb"))
            vprint("Saved TeamInfo Info",1)
            getTime(info.start_time)
            #info = cp.load(open("Temp/TeamInfo.info", "rb"))
            
            info.setIC(ic_path)
            
            cp.dump(info, open("Temp/IC.info","wb"))
            vprint("Saved IC Info",1)
            getTime(info.start_time)
            #info = cp.load( open("Temp/IC.info","rb"))
            
            info.setBenchmark(obocounts, benchmark_directory)
            
            cp.dump(info, open("Temp/Benchmark.info","wb"))
            vprint("Saved Benchmark Info",1)
            getTime(info.start_time)
            #info = cp.load( open("Temp/Benchmark.info","rb"))
            
            info.setPrediction()
            
            cp.dump(info, open("Temp/Prediction.info","wb"))
            vprint("Saved Prediction Info",1)
            getTime(info.start_time)
            #info = cp.load(open("Temp/Prediction.info","rb"))
            
            info.setInfo()
            
            cp.dump(info, open("Temp/Info.info","wb"))
            vprint("Saved Final Info",1)
            getTime(info.start_time)
            #info = cp.load(open("Temp/Info.info","rb"))
            
            
            vprint("Done prep work", 1)
            # For each mode
            for mode in ['partial', 'full']:
            #for mode in ['full']:
                # Threshold 0 should have the most proteins
                # If no coverage, there is no prediction. Skip!
                if(info.ProteinInPrediction[0.0] == 0 or info.Coverage is None):
                    vprint("No predicted proteins", 1)                    
                    continue
                
      
                
                info.setMode(mode)
                vprint("Mode: {}".format(mode), 1)
                getTime(info.start_time)
                # RUN METRICS
                # Set root path                
                summary_path = "{}/{}_{}_{}_Results".format(info.ResultPath, ontology, Type, mode)
                clear(summary_path)
                print("FMAX")
                F = FMAX(info)
                vwrite("FMAX: {}".format(F), summary_path, 1)
                
                print("WFMAX")
                WF = WFMAX(info)
                vwrite("WFMAX: {}".format(F), summary_path, 1)
                
                print("SMIN")
                S = SMIN(info)
                vwrite("SMIN: {}".format(F), summary_path, 1)
                
                print("NSMIN")
                NS = NSMIN(info)
                vwrite("NSMIN: {}".format(F), summary_path, 1)
                
                
                #TESTING AREA            
                vprint("Coverage: {}".format(info.Coverage), 1)
                vprint("ProteinInPrediction[0.00] : {}".format(info.ProteinInPrediction[0.00]), 1)
                vprint("ProteinInBenchmark : {}".format(info.ProteinInBenchmark), 1)
                
                vprint("Done", 1)
                getTime(info.start_time)
