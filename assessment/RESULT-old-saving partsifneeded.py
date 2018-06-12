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
        FMAXhandle  = open(p + "/{}_{}_{}_{}_{}{}_FMAX_results.txt".format(Ontology, self.taxon, Type, Mode, self.author, self.model), 'w')
        FMAXhandle.write('AUTHOR:    \t {} \n'.format(self.author))
        FMAXhandle.write('MODEL:   \t\t {} \n'.format(self.model))
        FMAXhandle.write('KEYWORDS:  \t {} \n'.format(self.keywords))
        FMAXhandle.write('SPECIES: \t\t {} \n'.format(self.taxon))
        
        FMAXhandle.write('{}\t\t   {}\t {}\t | {}         {}           {}    \n'.format('Ontology', 'Type', 'Mode', 'FMAX', 'Threshold', 'Coverage' ))
        FMAXhandle.write('{}\t\t\t {}\t {}\t | {:.2f}\t\t {:.2f}\t\t\t {:.2f}\n'.format( Ontology,  Type, Mode, self.FMAX, self.FMAXThreshold, coverage   ))
        
        FMAXhandle.write('{}\t | {}\t\t | {}\t\t | {}\n'.format('Threshold', 'PR', 'RC', 'F')) 
        FMAXhandle.write('------------------------------------------------------------------------------\n')
        index = 0
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            try:
                FMAXhandle.write('{:.2f}\t\t | {:.6f}\t | {:.6f}\t | {:.6f}\n'.format(threshold, self.PR[index], self.RC[index], self.F[index]))
            except IndexError:
                pass
            index += 1
        
        FMAXhandle.close()
        
        
        WFMAXhandle = open(p + "/%s_%s_%s_%s_%s%s_WFMAX_results.txt" % (Ontology, self.taxon, Type, Mode, self.author, self.model), 'w')
        WFMAXhandle.write('AUTHOR:    \t {} \n'.format(self.author))
        WFMAXhandle.write('MODEL:   \t\t {} \n'.format(self.model)) 
        WFMAXhandle.write('KEYWORDS:  \t {} \n'.format(self.keywords))  
        WFMAXhandle.write('SPECIES: \t\t {} \n'.format(self.taxon))
        WFMAXhandle.write('{}\t\t   {}\t {}\t | {}         {}           {}    \n'.format('Ontology', 'Type', 'Mode', 'WFMAX'   , 'Threshold'        , 'Coverage' ))
        WFMAXhandle.write('{}\t\t\t {}\t {}\t | {:.2f}\t\t {:.2f}\t\t\t {:.2f}\n'.format( Ontology ,  Type ,  Mode , self.WFMAX, self.WFMAXThreshold, coverage   ))
        
        WFMAXhandle.write('{}\t | {}\t\t | {}\t\t | {}\n'.format('Threshold', 'WPR', 'WRC', 'WF'))
        WFMAXhandle.write('------------------------------------------------------------------------------\n')
        index = 0
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            try:
                WFMAXhandle.write('{:.2f}\t\t | {:.6f}\t | {:.6f}\t | {:.6f}\n'.format(threshold  , self.WPR[index], self.WRC[index], self.WF[index]))
            except IndexError:
                pass
            index += 1
        
        WFMAXhandle.close()
        
        
        SMINhandle  = open(p + "/%s_%s_%s_%s_%s%s_SMIN_results.txt" % (Ontology, self.taxon, Type, Mode, self.author, self.model), 'w')
        SMINhandle.write('AUTHOR:    \t {} \n'.format(self.author))
        SMINhandle.write('MODEL:   \t\t {} \n'.format(self.model))
        SMINhandle.write('KEYWORDS:  \t {} \n'.format(self.keywords))  
        SMINhandle.write('SPECIES: \t\t {} \n'.format(self.taxon))
        SMINhandle.write('{}\t\t   {}\t {}\t | {}         {}           {}    \n'.format('Ontology', 'Type', 'Mode', 'SMIN'   , 'Threshold'       , 'Coverage' ))
        SMINhandle.write('{}\t\t\t {}\t {}\t | {:.2f}\t\t {:.2f}\t\t\t {:.2f}\n'.format( Ontology ,  Type ,  Mode , self.SMIN, self.SMINThreshold, coverage   ))
        
        SMINhandle.write('{}\t | {}\t\t | {}\t\t | {}\n'.format('Threshold', 'RU', 'MI', 'S'))
        SMINhandle.write('------------------------------------------------------------------------------\n')
        index = 0
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            try:
                SMINhandle.write('{:.2f}\t\t | {:.6f}\t | {:.6f}\t | {:.6f}\n'.format(threshold  , self.RU[index], self.MI[index], self.S[index]))
            except IndexError:
                pass
            index += 1
        
        SMINhandle.close()
        
        
        NSMINhandle = open(p + "/%s_%s_%s_%s_%s%s_NSMIN_results.txt" % (Ontology, self.taxon, Type, Mode, self.author, self.model), 'w')
        NSMINhandle.write('AUTHOR:    \t {} \n'.format(self.author))
        NSMINhandle.write('MODEL:   \t\t {} \n'.format(self.model))
        NSMINhandle.write('KEYWORDS:  \t {} \n'.format(self.keywords))  
        NSMINhandle.write('SPECIES: \t\t {} \n'.format(self.taxon))
        NSMINhandle.write('{}\t\t   {}\t {}\t | {}         {}           {}    \n'.format('Ontology', 'Type', 'Mode', 'NSMIN'   , 'Threshold'        , 'Coverage' ))
        NSMINhandle.write('{}\t\t\t {}\t {}\t | {:.2f}\t\t {:.2f}\t\t\t {:.2f}\n'.format( Ontology ,  Type ,  Mode , self.NSMIN, self.NSMINThreshold, coverage   ))
        
        NSMINhandle.write('{}\t | {}\t\t | {}\t\t | {}\n'.format('Threshold', 'NRU', 'NMI', 'NS'           ))
        NSMINhandle.write('------------------------------------------------------------------------------\n')
        index = 0
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            try:
                NSMINhandle.write('{:.2f}\t\t | {:.6f}\t | {:.6f}\t | {:.6f}\n'.format(threshold  , self.NRU[index], self.NMI[index], self.NS[index]))
            except IndexError:
                pass
            index += 1
        
        NSMINhandle.close()        
        
        