"""
A helper class for functions needed throughout the program 
"""

import os
import argparse
import sys
import yaml
import errno
import numpy
import time 



def printTime(start_time, message):
    '''
    Facilate Time broadcasts for optimizing
    '''
    
    if start_time == 0:
        start_time = time.time()
    
    print(message)
    print(time.time() - start_time)
    return start_time
    
    
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


def get_namespace_index(namespace):
    '''
    Convert namespace into indices
    
    Input:
    namespace : String    namespace
    
    Output:
    [0]       : Integer   corresponding value
    '''
    
    num = None
    if   namespace == 'BPO' or namespace == 'bpo':
        num = 0
    elif namespace == 'MFO' or namespace == 'mfo':
        num = 1
    elif namespace == 'CCO' or namespace == 'cco':
        num = 2
    elif namespace == 'HPO' or namespace == 'hpo':
        num = 2
    else:
        raise ValueError("name space not found, check prediction files")
        print(namespace)
    return num
    
     
def typeConverter(oldType):
    '''
    Description
    
    Input:
    oldType : String  old type name
    
    Output:
    [0]     : String  new type name
    '''
    
    if   oldType == 'type1':
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
    Base method for directory creation
    '''
    
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
            
            
'''
This is a nested group of methods that generate the entire directory structure 
for results and intermediate data
'''       
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

    
def read_config_MAIN():
    '''
    Read in the configuration file
    
    Output:
    [0] : String   OBO         file path
    [1] : String   IC          file path
    [2] : String   prediction  file path
    [3] : String   Benchmark   folder path
    [4] : String   Results     folder path
    [5] : String   verbose     Y or N   
    [6] :  
    -> If need more
    '''
    
    parser = argparse.ArgumentParser(description = 'Precision- Recall assessment for CAFA predictions.', )
    
    parser.add_argument('config_stream',type = extant_file, help = 'Configuration file')                
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
    
    return(obo_path, ic_path, prediction_path, benchmark_directory, results_directory)
    
    
def read_config_PLOT():
    '''
    NEEDS WORK
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
    title          = config_dict['title']
    smooth         = config_dict['smooth']
    if (smooth == 'Y'):
        s = True
    else:
        s = False
        
    methods = set()
    for i in xrange(Num_files):
        keyname = 'file'+str(i+1)
        methods.add(config_dict[keyname])
    return(results_folder, title, s, methods)   