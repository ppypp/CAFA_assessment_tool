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
        
        
        # FMAX
        self.FMAX           = 0.0
        self.PR             = []
        self.RC             = []
        self.F              = []
        self.FMAXThreshold  = 0.0
        
        # Weighted FMAX
        self.WFMAX          = 0.0
        self.WPR            = []
        self.WRC            = []
        self.WF             = []
        self.WFMAXThreshold = 0.0
        
        # SMIN
        self.SMIN           = 0.0
        self.RU             = []
        self.MI             = []
        self.S              = []
        self.SMINThreshold  = 0.0
        
        # Normalized SMIN        
        self.NSMIN          = 0.0
        self.NRU            = []
        self.NMI            = []
        self.NS             = []
        self.NSMINThreshold = 0.0
        
        # Raw data
        # Turth values (per protein?)
        self.TP             = 0 #= defaultdict()#Key : Protein, Value: List of Terms
        self.FN             = 0 #= defaultdict()#Key : Protein, Value: List of Terms
        self.FP             = 0 #= defaultdict()#Key : Protein, Value: List of Terms
        self.TN             = 0 #= defaultdict()#Key : Protein, Value: List of Terms
        self.coverage       = 0.0
        
        
        
        
        
    def update(self, evaluation, subval1, subval2, subval3, value, threshold):
        '''
        Stores the outpu of each metric in a single object
        
        Input:
        type      : String
        value     : Float
        subval1   : List[Float]
        subval2   : List[Float]
        threshold : Float
        '''
        if   evaluation == "FMAX":
            self.FMAX           = value
            self.PR             = subval1
            self.RC             = subval2
            self.F              = subval3
            self.FMAXThreshold  = threshold
        
        elif evaluation == "WFMAX":
            self.WFMAX          = value
            self.WPR            = subval1
            self.WRC            = subval2
            self.WF             = subval3
            self.WFMAXThreshold = threshold
            
        elif evaluation == "SMIN":
            self.SMIN           = value
            self.RU             = subval1
            self.MI             = subval2
            self.S              = subval3
            self.SMINThreshold  = threshold   
            
        elif evaluation == "NSMIN":
            self.NSMIN          = value
            self.NRU            = subval1
            self.NMI            = subval2
            self.NS             = subval3
            self.NSMINThreshold = threshold   
            
        else: 
            # Do Nothing
            return
            
    def info(self, author, model, keywords, taxon, predictionFile):
        '''
        Store the top level info
        '''        
        self.author         = author
        self.model          = model
        self.keywords       = keywords
        self.taxon          = taxon
        self.predictionFile = predictionFile
        
                
        
    def printOut(self):
        '''
        Print Data to console
        '''
        print('FMAX: %s\n' % self.FMAX)
        print('threshold giving FMAX: %s\n' % self.FMAXThreshold)
        print(self.F)
        print('WFMAX: %s\n' % self.WFMAX)
        print('threshold giving WFMAX: %s\n' % self.WFMAXThreshold)
        print(self.WF)
        print('SMIN: %s\n' % self.SMIN)
        print('threshold giving SMIN: %s\n' % self.SMINThreshold)
        print(self.S)
        print('NSMIN: %s\n' % self.NSMIN)
        print('threshold giving NSMIN: %s\n' % self.NSMINThreshold)
        print(self.NS)
        
        
    def writeOut(self, Ontology, Type, Mode, coverage):
        '''
        Write data to file
        
        We want a file per evaluation
        
        '''
        
        p = self.path
        FMAXhandle  = open(p + "/%s_%s_%s_%s_%s%s_FMAX_results.txt" % (Ontology, self.taxon, Type, Mode, self.author, self.model,), 'w')
        FMAXhandle.write('AUTHOR: %s   \n' % self.author)
        FMAXhandle.write('MODEL: %s    \n' % self.model) 
        FMAXhandle.write('KEYWORDS: %s \n' % self.keywords)  
        FMAXhandle.write('Species: %s  \n' % self.taxon)
        
        FMAXhandle.write('%s\t%s\t%s\t | %s\t%s\t%s\n' % ('Ontology', 'Type', 'Mode', 'FMAX'   , 'Threshold'       , 'Coverage' ))
        FMAXhandle.write('%s\t%s\t%s\t | %s\t%s\t%s\n' % ( Ontology ,  Type ,  Mode , self.FMAX, self.FMAXThreshold, coverage   ))
        
        FMAXhandle.write(    '%s\t | %s\t%s\t%s\t\n' % ('Threshold', 'PR'          , 'RC'          , 'F'          )) 
        index = 0
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            try:
                FMAXhandle.write('%s\t | %s\t%s\t%s\t\n' % (threshold  , self.PR[index], self.RC[index], self.F[index]))
            except IndexError:
                pass
            index += 1
        
        FMAXhandle.close()
        
        
        WFMAXhandle = open(p + "/%s_%s_%s_%s_%s%s_WFMAX_results.txt" % (Ontology, self.taxon, Type, Mode, self.author, self.model,), 'w')
        WFMAXhandle.write('AUTHOR: %s   \n' % self.author)
        WFMAXhandle.write('MODEL: %s    \n' % self.model) 
        WFMAXhandle.write('KEYWORDS: %s \n' % self.keywords)  
        WFMAXhandle.write('Species: %s  \n' % self.taxon)
        WFMAXhandle.write('%s\t%s\t%s\t | %s\t%s\t%s\n' % ('Ontology', 'Type', 'Mode', 'WFMAX'   , 'Threshold'        , 'Coverage' ))
        WFMAXhandle.write('%s\t%s\t%s\t | %s\t%s\t%s\n' % ( Ontology ,  Type ,  Mode , self.WFMAX, self.WFMAXThreshold, coverage   ))
        
        WFMAXhandle.write(    '%s\t | %s\t%s\t%s\t\n' % ('Threshold', 'WPR', 'WRC', 'WF'))
        index = 0
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            try:
                WFMAXhandle.write('%s\t | %s\t%s\t%s\t\n' % (threshold  , self.WPR[index], self.WRC[index], self.WF[index]))
            except IndexError:
                pass
            index += 1
        
        WFMAXhandle.close()
        
        
        SMINhandle  = open(p + "/%s_%s_%s_%s_%s%s_SMIN_results.txt" % (Ontology, self.taxon, Type, Mode, self.author, self.model,), 'w')
        SMINhandle.write('AUTHOR: %s   \n' % self.author)
        SMINhandle.write('MODEL: %s    \n' % self.model) 
        SMINhandle.write('KEYWORDS: %s \n' % self.keywords)  
        SMINhandle.write('Species: %s  \n' % self.taxon)
        SMINhandle.write('%s\t%s\t%s\t | %s\t%s\t%s\n' % ('Ontology', 'Type', 'Mode', 'SMIN'   , 'Threshold'       , 'Coverage' ))
        SMINhandle.write('%s\t%s\t%s\t | %s\t%s\t%s\n' % ( Ontology ,  Type ,  Mode , self.SMIN, self.SMINThreshold, coverage   ))
        
        SMINhandle.write(    '%s\t | %s\t%s\t%s\t\n' % ('Threshold', 'RU', 'MI', 'S'))
        index = 0
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            try:
                SMINhandle.write('%s\t | %s\t%s\t%s\t\n' % (threshold  , self.RU[index], self.MI[index], self.S[index]))
            except IndexError:
                pass
            index += 1
        
        SMINhandle.close()
        
        
        NSMINhandle = open(p + "/%s_%s_%s_%s_%s%s_NSMIN_results.txt" % (Ontology, self.taxon, Type, Mode, self.author, self.model,), 'w')
        NSMINhandle.write('AUTHOR: %s   \n' % self.author)
        NSMINhandle.write('MODEL: %s    \n' % self.model) 
        NSMINhandle.write('KEYWORDS: %s \n' % self.keywords)  
        NSMINhandle.write('Species: %s  \n' % self.taxon)
        NSMINhandle.write('%s\t%s\t%s\t | %s\t%s\t%s\n' % ('Ontology', 'Type', 'Mode', 'NSMIN'   , 'Threshold'        , 'Coverage' ))
        NSMINhandle.write('%s\t%s\t%s\t | %s\t%s\t%s\n' % ( Ontology ,  Type ,  Mode , self.NSMIN, self.NSMINThreshold, coverage   ))
        
        NSMINhandle.write(    '%s\t | %s\t%s\t%s\t\n' % ('Threshold', 'NRU', 'NMI', 'NS'           ))
        index = 0
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            try:
                NSMINhandle.write('%s\t | %s\t%s\t%s\t\n' % (threshold  , self.NRU[index], self.NMI[index], self.NS[index]))
            except IndexError:
                pass
            index += 1
        
        NSMINhandle.close()        
        
        