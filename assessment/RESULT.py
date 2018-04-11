# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 13:48:31 2018

@author: mcgerten
"""
import helper
import os
import numpy

##########################################REFACTOR#############################
# Its weird to have a class in a file by itself, should just reference the file 
# -> need to refactor references to this file
class result:
    ''' Stores results in a common format '''   
    
    
    def __init__(self, results_path):
        ''' 
        State all variables needed 
        
        Input:
        results_path : String       The directory to store results in
        '''
        self.path           = results_path         

        self.Opt            = 0.0
        self.Value1         = []
        self.Value2         = []
        self.Value3         = []
        self.OptThreshold   = 0.0
        
        self.coverage       = 0.0
        
        
        
        
        
    def update(self, subval1, subval2, subval3, value, threshold):
        '''
        Stores the outpu of each metric in a single object
        
        Input:
        type      : String
        value     : Float
        subval1   : List[Float]
        subval2   : List[Float]
        subval3   : List[Float]
        threshold : Float
        '''

        self.Opt            = float(value)
        self.Value1         = subval1
        self.Value2         = subval2
        self.Value3         = subval3
        self.OptThreshold   = float(threshold)
            
            
    def info(self, author, model, keywords, taxon, predictionFile):
        '''
        Store the top level info
        '''        
        self.author         = author
        self.model          = model
        self.keywords       = keywords
        self.taxon          = taxon
        self.predictionFile = predictionFile
        

    def writeOut(self, Ontology, Type, Mode, tool):
        '''
        
        '''
        
        # Set Ontology to correct formating
        if (Ontology == 'mfo'):
            O = 'MFO'
        elif (Ontology == 'bpo'):
            O = 'BPO'
        elif (Ontology == 'cco'):
            O ='CCO'
        else:
            O = Ontology
        # Set Type to correct format
        Type = helper.typeConverter(Type) 
        
        p = self.path
        handle  = open(p + "/{}_{}_{}_{}_{}{}_{}_results.txt".format(O, self.taxon, Type, Mode, self.author, self.model, tool), 'w')
        handle.write('!AUTHOR:      \t {}  \n'.format(self.author))
        handle.write('!MODEL:       \t {}  \n'.format(self.model))
        handle.write('!KEYWORDS:    \t {}  \n'.format(self.keywords))
        handle.write('!SPECIES:     \t {}  \n'.format(self.taxon))
        handle.write('!ONTOLOGY:    \t {}  \n'.format(O))
        handle.write('!TYPE:        \t {}  \n'.format(Type))
        handle.write('!MODE:        \t {}  \n'.format(Mode))
        handle.write('<OPTIMAL:\t {:.6f}   \n'.format(self.Opt))
        handle.write('<THRESHOLD:\t {:.2f} \n'.format(self.OptThreshold))
        handle.write('<COVERAGE:\t {:.2f}  \n'.format(self.coverage))

        
        handle.write('{}\t {}\t {}\t {}\n'.format('Threshold', 'PR', 'RC', 'F')) 
        index = 0
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            try:
                handle.write('>{:.2f}\t {:.6f}\t {:.6f}\t {:.6f}\n'.format(threshold, self.Value1[index], self.Value2[index], self.Value3[index]))
            except (IndexError, TypeError):
                pass
            index += 1
        
        handle.close()       
        
        
    