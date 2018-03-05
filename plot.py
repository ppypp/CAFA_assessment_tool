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


import helper


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
        Type   : String    {NK, LK} 
        mode   : String    {Full, Partial}
        method : String    Filename in full
        '''
        
        fields        = method.split('_')
        self.author   = fields[0]
        self.model    = int(fields[1])
        self.taxon    = fields[2]
        self.ontology = onto
        self.mode     = mode
        self.TYPE     = Type
        self.method   = method
        
        
    def read_values(self, evaluation, value1, value2, value3, opt, threshold, coverage):
        '''
        Read in the values 
        ############################################################### 
        SHould I just pass a filename and grab them directly?)
        ################################################################
        
        Input:
        evaluation : String        {FMAX, WFMAX, SMIN, NSMIN}
        value1     : List[Float]   Either PR or RU
        value2     : List[Float]   Either RC or MI
        value3     : List[Float]   Eval values for each val1/val2 pairing
        opt        : Float         Optimal Eval value
        threshold  : Float         Threshold opt occurs at
        coverage   : Float         Coverage 
        '''        
                
        if(evaluation == "FMAX"):
            self.FMAXPR       = value1
            self.FMAXRC       = value2
            self.FMAXV        = value3
            self.FMAXOPT      = opt
            self.FMAXTHR      = threshold
            self.coverage     = coverage
            
        elif(evaluation == "FMAX"):
            self.WFMAXPR      = value1
            self.WFMAXRC      = value2
            self.WFMAXV       = value3
            self.WFMAXOPT     = opt
            self.WFMAXTHR     = threshold
            self.coverage     = coverage
            
        elif(evaluation == "FMAX"):
            self.SMINRU       = value1
            self.SMINMI       = value2
            self.SMINV        = value3
            self.SMINOPT      = opt
            self.SMINTHR      = threshold
            self.coverage     = coverage
            
        elif(evaluation == "FMAX"):
            self.NSMINRU      = value1
            self.NSMINMI      = value2
            self.NSMINV       = value3
            self.NSMINOPT     = opt
            self.NSMINTHR     = threshold
            self.coverage     = coverage
            
        else:
            # Not a valid evaluation 
            pass
        
        
###############################################END OF RESULT ############################        
def getResults(title, method, onto, Type, mode, evaluation, results_folder):
    '''
    
    '''
    
    # Type need to be 'NK' and 'LK'
    filename = os.path.join(results_folder,'%s_%s.txt' % (method, evaluation))
    r = result()
    r.read_info(onto, Type, mode, method)
    with open(filename,'r') as f:
        read = False
        for line in f:
            if line.startswith('>'):
                fields = line.strip()[1:].split('\t')
                if fields[0] == onto and fields[1] == Type and fields[2] == mode:
                    read = True
                    #print(line)
                    v= f.next().split('|')[1].strip().split()
                   # print(pr)
                    v2 = f.next().strip().split()
                    #print(rc)
                    value1    = []
                    value2    = []
                    value3    = []
                    opt       = 0.0
                    threshold = 0.0
                    coverage  = 0.0
                    r.read_values(evaluation, value1, value2, value3, opt, threshold, coverage)
                    
        if not read:
            print('result not found for %s %s %s %s' % (method, onto, Type, mode))
    return(r)
    


def curveSmooth(v1, v2):
    '''
    Curve smoothing of a given curve
    
    Input:
    result     : Object
    
    Ouput:
    [0]        : List[float,float]   [Value 1, Value 2]
    '''
    
    #This function removes a pair if there exists another pair that's greater in both values
    #values should both be lists of the same length
    val1 = []
    val2 = []
     
        
    for i in range(len(v1)):
        remove = False
        for j in range(i):
            if v1[i] < v1[j] and v2[i] < v2[j]:
                remove = True
                break
        if not remove:
            val1.append(v1[i])
            val2.append(v2[i])
    return([val1, val2])
    
    
def plotMultiple(title, evaluation, results, smooth):
    '''
    Plots the graph
    
    Input:
    title   : String
    results : List[Objects]
    Smooth  : Boolean
    '''
    
    fontP  = FontProperties()
    fontP.set_size('small')
    num    = len(results)
    pal    = sns.color_palette("Paired", num)
    colors = pal.as_hex()
        
    
     
        
    for j, result in enumerate(results):
        linetype = '-'
        method   = result.method
        coverage = result.coverage
        # Grab values corresponding to evaluation type
        if   (evaluation == "FMAX"):
            value1      = results.FMAXPR
            value2      = results.FMAXRC
            opt         = results.FMAXOPT
            thres       = results.FMAXTHR
        elif (evaluation == "WFMAX"):
            value1      = results.WFMAXPR
            value2      = results.WFMAXRC
            opt         = results.WFMAXOPT
            thres       = results.WFMAXTHR
        elif (evaluation == "SMIN"):
            value1      = results.SMINRU
            value2      = results.SMINMI
            opt         = results.SMINOPT
            thres       = results.SMINTHR
        elif (evaluation == "NSMIN"): 
            value1      = results.NSMINRU
            value2      = results.NSMINMI
            opt         = results.NSMINOPT
            thres       = results.NSMINTHR
        else:
            # If not valid type, return empty Lists
            value1 = []
            value2 = []
        if smooth:
            ax = plt.subplot()
            v1 = curveSmooth(result, evaluation)[0][1:]
            v2 = curveSmooth(result, evaluation)[1][1:]
            ax.plot(v2, v1, linetype, color = colors[j], label = method + ':\nF=%s C=%s'%(opt, coverage)) 
            ax.plot(value2[int(thres*100)], value1[int(thres*100)], 'o', color = colors[j])
        else:
            ax = plt.subplot()
            ax.plot(value2, value1, linetype, color = colors[j], label = method + ':\nF=%s C=%s'%(opt, coverage))
            ax.plot(value2[int(thres*100)], value1[int(thres*100)], 'o', color = colors[j])
    plt.axis([0, 1, 0, 1 ])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    plt.yticks(numpy.arange(0, 1, 0.1))
    plt.xlabel('Val1')
    plt.ylabel('Val2')
    ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
    plt.title(title)
    figurename = os.path.join('./plots/',title)       
    plt.savefig(figurename, dpi = 200)
    plt.close()
        

def check_existence(results_folder, methods):
    '''
    
    '''
    
    # Assume exists until shown otherwise
    re = True
    if os.path.isdir(results_folder):
        for method in methods:
            file_result = os.path.join(results_folder,'%s_results.txt' % method)
            if not os.path.isfile(file_result):
                sys.stderr.write('file %s not found\n' % file_result)
                re = False
                break
    else:
        sys.stderr.write('directory %s not found\n' % results_folder)
        re = False
    return(re)
    
     
if __name__=='__main__':
    '''
    
    '''
    
    # Get config details
    results_folder, title, smooth, methods = helper.read_config_PLOT()
    if check_existence(results_folder, methods):
        rawdata_folder = results_folder + '/rawdata/'
        for evaluation in ['FMAX', 'WFMAX', 'SMIN', 'NSMIN']:
            for onto in ['BPO','CCO','MFO']:
                for Type in ['LK','NK']:
                    for mode in ['partial','full']:
                        specific_title = '%s_%s_%s_%s_%s.png' % (title, onto, Type, mode, evaluation)
                        print('\nPlotting %s\n' % title)
                        helper.mkdir_p('./plots')
                        result_list = []
                        for method in methods:
                            r = getResults(title, method, onto, Type, mode, evaluation, results_folder)
                            result_list.append(r)
                        plotMultiple(specific_title, result_list, smooth)
    else:
        sys.exit()
 
