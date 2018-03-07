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
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import sys
from matplotlib.font_manager import FontProperties


import helper


class result:    
    
    def __init__(self):
        '''
        Intilize result object
        '''
        
        self.evaluation   = ''
        self.Value0       = []
        self.Value1       = []
        self.Value2       = []
        self.Value3       = []
        self.Opt          = float
        self.Threshold    = float        
                
        self.author       = ''
        self.model        = ''
        self.keywords     = ''
        self.taxon        = ''
        self.ontology     = ''
        self.mode         = ''
        self.Type         = ''
        
        self.Coverage     = float
        
        
###############################################END OF RESULT ############################        
def getResults(method, onto, taxon, Type, mode, evaluation, results_folder):
    '''
    
    '''
    
    r = result()
    # Intialize Header    
    fields     = method.split('_')
    r.author   = fields[0]
    r.model    = int(fields[1])
    r.taxon    = taxon
    r.ontology = onto
    r.mode     = mode
    r.Type     = Type
    r.method   = method
    am = r.author + fields[1]
    
        
    
    filename = os.path.join(results_folder,'{}_{}_{}_{}_{}_{}_results.txt'.format(onto, taxon, Type, mode, am, evaluation))
    
    with open(filename,'r') as f:
        read      = False
        tvalue    = []
        value1    = []
        value2    = []
        value3    = []
        opt       = 0.0
        threshold = 0.0
        coverage  = 0.0
        
        for line in f:
            if(line.startswith('<')):
                read = True
                fields = line.strip().split('\t')
                if(fields[0]   == '<OPTIMAL:'):
                    opt        = fields[1] 
                elif(fields[0] == '<THRESHOLD:'):
                    threshold  = fields[1]
                elif(fields[0] == '<COVERAGE:'):
                    coverage   = fields[1]
            elif(line.startswith('>')):
                fields = line.strip()[1:].split('\t')
                tvalue.append(fields[0])
                value1.append(fields[1])
                value2.append(fields[2])
                value3.append(fields[3])
        if not read:
            print('result not found for %s %s %s %s' % (method, onto, Type, mode))
    # Set result values
    r.Value0    = tvalue # Thresholds
    r.Value1    = value1 # PR or RU
    r.Value2    = value2 # RC or MI
    r.Value3    = value3 # F or S
    r.Opt       = opt
    r.Threshold = threshold
    r.Coverage  = coverage
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
    
    
def plotMultiple(title, results, smooth):
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
        method      = result.method
        coverage    = result.Coverage
        value1      = result.Value3
        value2      = result.Value0
        opt         = result.Opt
        thres       = result.Threshold
        
        if smooth: # BROKEN 
            ax = plt.subplot()
            v1 = curveSmooth(result)[0][1:]
            v2 = curveSmooth(result)[1][1:]
            ax.plot(v2, v1, linetype, color = colors[j], label = method + ':\nF=%s C=%s'%(opt, coverage)) 
            ax.plot(value2[int(thres*100)], value1[int(thres*100)], 'o', color = colors[j])
        else:
            ax = plt.subplot()
            ax.plot(value2, value1, color = colors[j], label = method + ':\nOpt ={}\nCov ={}'.format(opt, coverage))
            #ax.plot(value2[int(thres*100)], value1[int(thres*100)], 'o',     color = colors[j])
    plt.autoscale(True,'both')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    # Use evaluation to label axises
    plt.xlabel('Threshold')
    plt.ylabel('Value')
    ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5))
    plt.title(title)
    figurename = os.path.join('./plots/', title)       
    plt.savefig(figurename, dpi = 300)
    plt.close()
        

     
if __name__=='__main__':
    '''
    
    '''
    
    # Get config details
    results_folder, title, smooth, methods = helper.read_config_PLOT()
    rawdata_folder = results_folder + '/rawdata/'
    for evaluation in ['FMAX', 'WFMAX', 'SMIN', 'NSMIN']:
        for onto in ['BPO','CCO','MFO']:
            for Type in ['LK','NK']:
                for mode in ['partial','full']:
                    specific_title = '%s_%s_%s_%s_%s.png' % (onto, title, Type, mode, evaluation)
                    print('\nPlotting %s\n' % title)
                    helper.mkdir_p('./plots')
                    result_list = []
                    for method in methods:
                        fields     = method.split('_')
                        taxon    = fields[2]
                        r = getResults(method, onto, taxon, Type, mode, evaluation, results_folder)
                        result_list.append(r)
                    plotMultiple(specific_title, result_list, smooth)
    else:
        sys.exit()
 
