#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from assessment.assessmentTools import Info, read_benchmark
from assessment.GOPrediction import GOPrediction
import os
import sys
import errno    
import gc
import yaml
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
 

def typeConverter(oldType):
    '''
    Description
    
    Input:
    oldType : String  old type name
    
    Output:
    [0]     : String  new type name
    '''
    
    if oldType=='type1':
        newType = 'NK'
    elif oldType == 'type2':
        newType = 'LK'
    elif oldType == 'all':
        newType = 'All'
    return(newType)
    

def extant_file(x):
    '''
    Description - extant?
    
    Input:
    x   : 
    
    Output:
    [0] :
    '''
    
    if not os.path.isfile(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    else:
        return(open(x,'r'))
        
        
def mkdir_p(path):
    '''
    Make a new directory checking for OS errors    
    
    Input:
    path : String    The filepath for the new directory
    '''
    
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            # If another error
            raise
            
            
def read_config():
    '''
    Read in the configuration file
    
    Output:
    [0] : String   OBO         file path
    [1] : String   Benchmark   folder path
    [2] : String   Results     folder path
    [3] : String   prediction  file path
    [4] : String   IC          file path
    [5] : String   verbose     Y or N   
    [6] : String   
    -> If need more
    '''
    
    parser = argparse.ArgumentParser(description='Precision- Recall assessment for CAFA predictions.', )
    
    parser.add_argument('config_stream',type=extant_file, help='Configuration file')
    # CAFA3 raw submission filename formats are listed here:https://www.synapse.org/#!Synapse:syn5840147/wiki/402192
    # example filename format: Doegroup_1_9606.txt/Doegroup_2_hpo.txt
    # If prediction file is already split by ontology it should follow Doegroup_1_9606_BPO.txt(or _MFO, _CCO)                  
    args = parser.parse_args()
    # Load config file to dictionary
    try:
        config_dict = yaml.load(args.config_stream)['assessment']
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()
    # Store config into itermediate variables    
    obo_path         = config_dict['obo_path']
    benchmark_path   = config_dict['benchmark_path']
    results_path     = config_dict['results_path']
    f                = config_dict['file']
    ic_path          = config_dict['ic_path']
    verbose          = config_dict['verbose']
    
    return(obo_path, benchmark_path, results_path, f, ic_path, verbose)


if __name__=='__main__':
    '''
    Main function that takes a predicition and returns calculated values
    '''
    # Read Config
    obo_path, benchmarkFolder, resultsFolder, f , ic_path, verbose = read_config()
    # Setup global variables
    helper.init(verbose)
    # Setup workspace
    mkdir_p(resultsFolder)
    mkdir_p(resultsFolder+'/rawdata/')
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
          
    resulthandle= open(resultsFolder + "/%s_results.txt" % (os.path.basename(f).split('.')[0]),'w')
    prhandle = open(resultsFolder + "/%s_prrc.txt" % (os.path.basename(f).split('.')[0]),'w')
    
    resulthandle.write('AUTHOR:%s\n' % author)
    resulthandle.write('MODEL: %s\n' % model) 
    resulthandle.write('KEYWORDS: %s\n' % keywords)  
    resulthandle.write('Species:%s\n' % taxon)
    resulthandle.write('%s\t%s\t%s\t | %s\t%s\t%s\n' % ('Ontology','Type','Mode','Fmax','Threshold','Coverage'))
    
    # Grab calculated IC
    ic = cp.load(open(ic_path,"rb"))
    # Make a result object for storing output
    r = result()
    
    for ontology in ['bpo','cco','mfo']:
        path = os.path.splitext(prediction_path.name)[0] + '_'+ontology.upper() + '.txt'
        print('ontology: %s\n' % ontology)
        for Type in ['type1','type2']:
            print('benchmark type:%s\n' % typeConverter(Type))
            benchmark, obocountDict = read_benchmark(ontology, taxon_name_converter(taxon), Type, benchmarkFolder, obo_path)
            if benchmark == None:
                sys.stderr.write('No benchmark is available for the input species and type')
            # Create information object that passing necessary information to subroutines    
            i = Info(benchmark, path, obocountDict[ontology], ic)
            # Check for success
            if i.exist:
                for mode in ['partial', 'full']:
                    print('mode:%s\n' % mode)
                    fm = i.check("FMAX", mode)
                    r.update("FMAX", fm[0], fm[1], fm[2], fm[3])
                    wfm = i.check("WFMAX", mode)
                    r.update("WFMAX", wfm[0], wfm[1], wfm[2], wfm[3])
                    sm = i.check("SMIN", mode)
                    r.update("SMIN", sm[0], sm[1], sm[2], sm[3])
                    nsm = i.check("NSMIN", mode)
                    r.update("NSMIN", nsm[0], nsm[1], nsm[2], nsm[3])
                    
                    coverage = i.coverage()
                    #fm.append(os.path.splitext(os.path.basename(pred_path.name))[0])
                    r.printOut()
                    print('coverage: %s\n' % coverage)
                    r.writeOut() # Pass in the handle -> Is this acceptable
                    
            del i
            gc.collect()
    resulthandle.close()
    prhandle.close()
