# -*- coding: utf-8 -*-
"""
A helper class for functions needed throughout the program 
"""


import os
import argparse
import sys
import yaml
import errno
import numpy

def init(v_choice):
    '''
    Setup the global variables needed

    Input:
    v_choice : String    whether to be verbose
    '''
    
    global v
    
    if(v_choice == "Y"):
        v = 1
    elif(v_choice == "N"):
        v = 0
    else:
        # Throw an error?
        pass


# Could make this a tiered system -> differing numbers represent different verbosity
def vprint(s):
    '''
    Print only if Verbose is chosen
    
    Input:
    s : The string to print
    '''
    if v == 1 :
        print(s)
     
     
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
    
    '''
    
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
            
def mkdir_results(path):
    mkdir_p(path + '/FMAX/')
    mkdir_ontology(path + '/FMAX/')
    mkdir_p(path + '/WFMAX/')
    mkdir_ontology(path + '/WFMAX/')
    mkdir_p(path + '/SMIN/')
    mkdir_ontology(path + '/SMIN/')
    mkdir_p(path + '/NSMIN/')
    mkdir_ontology(path + '/NSMIN/')
    mkdir_p(path + '/AUC/')
    mkdir_ontology(path + '/AUC/')
    
                   
def mkdir_ontology(path):
    mkdir_p(path + '/bpo/')
    mkdir_type(path + '/bpo/')
    mkdir_p(path + '/cco/')
    mkdir_type(path + '/cco/')
    mkdir_p(path + '/mfo/')
    mkdir_type(path + '/mfo/')
    mkdir_p(path + '/hpo/')
    mkdir_type(path + '/hpo/')

def mkdir_type(path):
    mkdir_p(path + '/type1/')
    mkdir_mode(path + '/type1/')
    mkdir_p(path + '/type2/')
    mkdir_mode(path + '/type2/')
    
def mkdir_mode(path):   
    mkdir_p(path + '/full/')
    mkdir_threshold(path + '/full/')
    mkdir_p(path + '/partial/')
    mkdir_threshold(path + '/partial/')
    
def mkdir_threshold(path):
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
        threshold = numpy.around(threshold, decimals = 2)
        mkdir_p('{}/{}/'.format(path, threshold))
    
    
          
def read_config_IC():
    '''
    Read in the configuration file for IC Tool
    
    Output:
    [0] : String   OBO         file path
    [1] : String   GAF         file path
    '''
    
    parser = argparse.ArgumentParser(description='Precision- Recall assessment for CAFA predictions.', )
    
    parser.add_argument('config_stream', type = extant_file, help = 'Configuration file')            
    args = parser.parse_args()
    # Load config file to dictionary
    try:
        config_dict = yaml.load(args.config_stream)['assessment']
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()
    # Store config into itermediate variables    
    obo_path = config_dict['obo_path']
    gaf_path = config_dict['gaf_path']
     
    return(obo_path, gaf_path)
    
    
def read_config_MAIN():
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
    
    parser = argparse.ArgumentParser(description = 'Precision- Recall assessment for CAFA predictions.', )
    
    parser.add_argument('config_stream',type = extant_file, help = 'Configuration file')
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
    obo_path             = config_dict['obo_path']
    ic_path              = config_dict['ic_path']
    prediction_path      = config_dict['prediction_path']
    benchmark_directory  = config_dict['benchmark_path']
    results_directory    = config_dict['results_path']
    verbose              = config_dict['verbose']
    
    return(obo_path, ic_path, prediction_path, benchmark_directory, results_directory, verbose)
    
    
def read_config_PLOT():
    '''
    
    '''
    
    parser = argparse.ArgumentParser(description='Precision-Recall Curves plot', )
    

    parser.add_argument('config_stream',type=extant_file, help='Configuration file')
    #CAFA3 raw submission filename formats are listed here:https://www.synapse.org/#!Synapse:syn5840147/wiki/402192
    #example filename format: Doegroup_1_9606.txt/Doegroup_2_hpo.txt
    #If prediction file is already split by ontology it should follow Doegroup_1_9606_BPO.txt(or _MFO, _CCO)                  
    args = parser.parse_args()
    try:
        config_dict = yaml.load(args.config_stream)['plot']
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()
    Num_files = len(config_dict)-3
    results_folder = config_dict['results']
    title = config_dict['title']
    smooth = config_dict['smooth']
    if (smooth == 'Y'):
        s = True
    else:
        s = False
        
    methods = set()
    for i in xrange(Num_files):
        keyname = 'file'+str(i+1)
        methods.add(config_dict[keyname])
    return(results_folder, title, s, methods)   