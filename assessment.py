#!/usr/bin/env python
# -*- coding: utf-8 -*-

from assessment.assessmentTools import Info, read_benchmark
from assessment.GOPrediction import GOPrediction
import os
import sys
import gc
import pickle as cp
from assessment.RESULT import result

import helper


def get_namespace_index(namespace):
    '''
    Convert namespace into indices
    
    Input:
    namespace : String    namespace
    
    Output:
    [0]       : Integer   corresponding value
    '''
    
    num = None
    if namespace=='BPO' or namespace=='bpo':
        num = 0
    elif namespace=='MFO' or namespace=='mfo':
        num = 1
    elif namespace=='CCO' or namespace=='cco':
        num =2
    else:
        raise ValueError("name space not found, check prediction files")
        print(namespace)
    return num


def taxon_name_converter(taxonID):
    '''
    Convert from taxonomy ID to name (i.e. from 9606 to HUMAN）
    
    Input:
    taxonID : String     representing ID
    
    Output:
    [0]     : String     the correspnding taxon name
    '''

    taxonTable = {
    '10116':'RAT',    '9606':'HUMAN',   '3702':'ARATH',   '7955':'DANRE',
    '44689':'DICDI',  '7227':'DROME',   '83333':'ECOLI',  '10090':'MOUSE',
    '208963':'PSEAE', '237561':'CANAX', '559292':'YEAST', '284812':'SCHPO',
    '8355':'XENLA',   '224308':'BACSU', '99287':'SALTY',  '243232':'METJA',
    '321314':'SALCH', '160488':'PSEPK', '223283':'PSESM', '85962':'HELPY',
    '243273':'MYCGE', '170187':'STRPN', '273057':'SULSO', 
    'all':'all', 'prokarya':'prokarya', 'eukarya':'eukarya'}
    
    return taxonTable[taxonID]    


if __name__=='__main__':
    '''
    Main function that takes a predicition and returns calculated values
    '''
    # Read Config
    obo_path, benchmarkFolder, resultsFolder, f , ic_path, verbose = helper.read_config_MAIN()
    # Setup global variables
    helper.init(verbose)
    # Setup workspace
    helper.mkdir_p(resultsFolder)
    helper.mkdir_p(resultsFolder+'/rawdata/')
    print('\nEvaluating %s.\n' % f)
    # Get predictions
    all_prediction = GOPrediction()
    prediction_path = open(f, 'r')
    all_prediction.read_and_split_and_write(obo_path, prediction_path)
    info = [all_prediction.author, all_prediction.model, all_prediction.keywords, all_prediction.taxon]
    # Clear memory
    del all_prediction
    gc.collect()
    # Store values
    author = info[0]
    model = info[1]
    keywords = info[2][0] 
    taxon = info[3]
    print('AUTHOR: %s\n' % author)
    print('MODEL: %s\n' % model)
    print('KEYWORDS: %s\n' % keywords)
    print('Species:%s\n' % taxon)
    
    # Grab calculated IC
    ic = cp.load(open(ic_path,"rb"))
    # Make a result object for storing output
    r = result(resultsFolder)
    r.info(author, model, keywords, taxon, f)
    
    for ontology in ['bpo','cco','mfo']:
        path = os.path.splitext(prediction_path.name)[0] + '_'+ontology.upper() + '.txt'
        print('ontology: %s\n' % ontology)
        for Type in ['type1','type2']:
            print('benchmark type:%s\n' % helper.typeConverter(Type))
            benchmark, obocountDict = read_benchmark(ontology, taxon_name_converter(taxon), Type, benchmarkFolder, obo_path)
            if benchmark == None:
                sys.stderr.write('No benchmark is available for the input species and type')
            # Create information object that passing necessary information to subroutines    
            i = Info(benchmark, path, obocountDict[ontology], ic)
            # Check for success
            if i.exist:
                for mode in ['partial', 'full']:
                    
                    fm = i.check("FMAX", mode)
                    r.update("FMAX", fm[0], fm[1], fm[2], fm[3], fm[4])
                    wfm = i.check("WFMAX", mode)
                    r.update("WFMAX", wfm[0], wfm[1], wfm[2], wfm[3], wfm[4])
                    sm = i.check("SMIN", mode)
                    r.update("SMIN", sm[0], sm[1], sm[2], sm[3], sm[4])
                    nsm = i.check("NSMIN", mode)
                    r.update("NSMIN", nsm[0], nsm[1], nsm[2], nsm[3], nsm[4])
                    coverage = i.coverage()
                    # Write to file
                    r.writeOut(ontology, Type, mode, coverage) 
                    
            del i
            gc.collect()
