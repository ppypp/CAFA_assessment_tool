#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 11:32:54 2017
This script provides top ten labs from the rankings and output their respective prediction files
@author: nzhou
"""

import os
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy
import sys
from matplotlib.font_manager import FontProperties
import yaml
import argparse
import errno
import pickle as cp



def typeConverter(oldType):
    '''
    Converts from outdated types to current ones.
    
    Input:
    oldtype : String    old version of type
    
    Output:
    [0]     : String    newtype
    '''
    if oldType=='type1':
        newType = 'NK'
    elif oldType == 'type2':
        newType = 'LK'
    elif oldType == 'all':
        newType = 'All'
    return(newType)


class result:    
    
    def __init__(self):
        '''
        Intilize result object
        '''
        
        self.type         = ''
        #type include: precision/recall, weighted pr, ru/mi, weighted ru/mi
        self.FMAXPR       = []
        self.FMAXRC       = []
        self.FMAXOPT      = float
        self.FMAXTHR      = float
        
        self.WFMAXPR      = []
        self.WFMAXRC      = []
        self.WFMAXOPT     = float
        self.WFMAXTHR     = float
        
        self.SMINRU       = []
        self.SMINMI       = []
        self.SMINOPT      = float
        self.SMINTHR      = float
                
        self.NSMINRU      = []
        self.NSMINMI      = []
        self.NSMINOPT     = float
        self.NSMINTHR     = float
        
        self.author       = ''
        self.model        = ''
        self.keywords     = ''
        self.taxon        = ''
        self.ontology     = ''
        self.mode         = ''
        self.TYPE         = ''
        
        self.coverage     = float
        #NK or LK
        
        
    def read_info(self, onto, Type, mode, method):
        '''
        Store data
        
        Input:
        onto   : String    Ontology
        Type   : 
        mode   : String    {Full, Partial}
        method : 
        '''
        
        fields        = method.split('_')
        self.author   = fields[0]
        self.model    = int(fields[1])
        self.taxon    = fields[2]
        self.ontology = onto
        self.mode     = mode
        self.TYPE     = Type
        self.method   = method
        
        
    def read_values(self, evaluation, value1, value2, opt, threshold):
        '''
        Read in the values 
        ############################################################### 
        SHould I just pass a filename and grab them directly?)
        ################################################################
        
        Input:
        evaluation : String        {FMAX, WFMAX, SMIN, NSMIN}
        value1     : List[Float]
        value2     : List[Float]
        '''        
                
        if(evaluation == "FMAX"):
            self.FMAXPR       = value1
            self.FMAXRC       = value2
            self.FMAXOPT      = opt
            self.FMAXTHR      = threshold
            
        elif(evaluation == "FMAX"):
            self.WFMAXPR      = value1
            self.WFMAXRC      = value2
            self.WFMAXOPT     = opt
            self.WFMAXTHR     = threshold
            
        elif(evaluation == "FMAX"):
            self.SMINRU       = value1
            self.SMINMI       = value2
            self.SMINOPT      = opt
            self.SMINTHR      = threshold
            
        elif(evaluation == "FMAX"):
            self.NSMINRU      = value1
            self.NSMINMI      = value2
            self.NSMINOPT     = opt
            self.NSMINTHR     = threshold
            
        else:
            # Not a valid evaluation 
            pass
    
            
    def getCoverage(self,onto,Type,mode,method,results_folder):
        '''
        Needed?
        '''
        
        result_file = os.path.join(results_folder,'%s_results.txt' % method)
        with open(result_file) as f:
            for line in f:
                if line.startswith(onto):
                    if line.split('|')[0].split()[1]==Type and line.split('|')[0].split()[2]==mode:
                        coverage = line.split('|')[1].split()[2]
                        self.coverage=str(numpy.around(float(coverage),decimals=2))
                        break            
        
        
###############################################END OF RESULT ############################        
def getprrc(onto,Type,mode,prrcfolder,method, results_folder):
    '''
    
    '''
    
    #Type need to be 'NK' and 'LK'
    filename = os.path.join(prrcfolder,'%s_prrc.txt' % (method))
    r = result()
    r.read_info(onto,Type,mode,method)
    with open(filename,'r') as f:
        read = False
        for line in f:
            if line.startswith('>'):
                fields = line.strip()[1:].split('\t')
                if fields[0]==onto and fields[1]==Type and fields[2]==mode:
                    read=True
                    #print(line)
                    pr = f.next().split('|')[1].strip().split()
                   # print(pr)
                    rc = f.next().strip().split()
                    #print(rc)
                    r.read_prrc(pr,rc)
                    r.calculate_fmax()
                    if r.check_fmax(onto,Type,mode,method, results_folder):
                        r.getCoverage(onto,Type,mode,method, results_folder)                
                    else:
                        print("check results\n")
                    break
                #both pr and rc are lists
        if not read:
            print('result not found for %s %s %s %s' % (method, onto, Type, mode))
    return(r)
    


def curveSmooth(result, t):
    '''
    Curve smoothing of a given curve
    
    Input:
    result : Object
    t      : String              Type of ansyslis
    Ouput:
    [0]    : List[float,float]   [Precision, Recall]
    '''
    
    #This function removes a p-r pair if there exists another p-r pair that's greater in both precision and recall
    #precision and recall should both be lists of the same length
    val1 = []
    val2 = []
    val1name = ""
    val2name = ""
    
    if (t == "FMAX"):
        val1name = "FMAXPR"
        val2name = "FMAXRC"
    elif (t == "WFMAX"):
        val1name = "WFMAXPR"
        val2name = "WFMAXRC"
    elif (t == "SMIN"):
        val1name = "SMINRU"
        val2name = "SMINMI"
    elif(t == "NSMIN"): 
        val1name = "NSMINRU"
        val2name = "NSMINMI"
    else:
        # If not valid type, return empty Lists
        return([[],[]])        
        
    for i in range(len(result.val1name)):
        remove = False
        for j in range(i):
            if result.val1name[i]<result.val1name[j] and result.val2name[i]<result.val2name[j]:
                remove = True
                break
        if not remove:
            val1.append(result.val1name[i])
            val2.append(result.val2name[i])
    return([val1,val2])
    
    
def plotMultiple(title, listofResults, smooth):
    '''
    Plots the graph
    
    Input:
    title   : String
    results : List[]
    Smooth  : Boolean
    '''
    
    fontP = FontProperties()
    fontP.set_size('small')
    num = len(listofResults)
    pal=sns.color_palette("Paired", num)
    colors=pal.as_hex()
    for j,i in enumerate(listofResults):
        linetype = '-'
        if smooth=='Y':
            ax = plt.subplot()
            precision = curveSmooth(i)[0][1:]
            recall = curveSmooth(i)[1][1:]
            ax.plot(recall,precision,linetype,color=colors[j],label=i.method+':\nF=%s C=%s'%(i.opt,i.coverage)) 
            ax.plot(i.recall[int(i.thres*100)],i.precision[int(i.thres*100)],'o',color=colors[j])
        elif smooth=='N':
            ax = plt.subplot()
            ax.plot(i.recall,i.precision,linetype,color=colors[j],label=i.method+':\nF=%s C=%s'%(i.opt,i.coverage))
            ax.plot(i.recall[int(i.thres*100)],i.precision[int(i.thres*100)],'o',color=colors[j])
    plt.axis([0,1,0,1])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    plt.yticks(numpy.arange(0,1,0.1))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title(title)
    figurename = os.path.join('./plots/',title)       
    plt.savefig(figurename,dpi=200)
    plt.close()
        
        
        
def extant_file(x):
    '''
    
    '''
    
    if not os.path.isfile(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    else:
        return(open(x,'r')) 
        
        
def read_config():
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
    methods = set()
    for i in xrange(Num_files):
        keyname = 'file'+str(i+1)
        methods.add(config_dict[keyname])
    return(results_folder, title,smooth, methods)   


def check_existence(results_folder, methods):
    '''
    
    '''
    
    re = True
    if os.path.isdir(results_folder):
        for method in methods:
            file_result = os.path.join(results_folder,'%s_results.txt' % method)
            if not os.path.isfile(file_result):
                sys.stderr.write('file %s not found\n' % file_result)
                re = False
                break
        prrc_folder = results_folder + '/pr_rc/'
        if os.path.isdir(prrc_folder):
            #check if files exist
            for method in methods:
                file_prrc = os.path.join(prrc_folder,'%s_prrc.txt' % method)
                if not os.path.isfile(file_prrc):
                    re = False
                    sys.stderr.write('file %s not found\n' % file_prrc)
                    break
        else:
            sys.stderr.write('directory %s not found\n' % prrc_folder)
            re = False
    else:
        sys.stderr.write('directory %s not found\n' % results_folder)
        re = False
    return(re)
    
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
     
if __name__=='__main__':
    '''
    
    '''
    
    # Get config details
    results_folder, title, smooth, methods = read_config()
    if check_existence(results_folder, methods):
        rawdata_folder = results_folder + '/rawdata/'
        for onto in ['bpo','cco','mfo']:
            for Type in ['LK','NK']:
                for mode in ['partial','full']:
                    specific_title = '%s_%s_%s_%s_fmax.png' % (title, onto, Type, mode)
                    print('\nPlotting %s\n' % title)
                    mkdir_p('./plots')
                    result_list = []
                    for method in methods:
                        res = getprrc(onto, Type, mode, prrcfolder, method, results_folder)
                        result_list.append(res)
                    plotMultiple(specific_title,result_list,smooth)
    else:
        sys.exit()
 
