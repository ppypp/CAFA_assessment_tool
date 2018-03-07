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
        
        if (Ontology == 'mfo'):
            O = 'MFO'
        elif (Ontology == 'bpo'):
            O = 'BPO'
        elif (Ontology == 'cco'):
            O ='CCO'
        else:
            O = Ontology
        Type = helper.typeConverter(Type) 
        
        p = self.path
        FMAXhandle  = open(p + "/{}_{}_{}_{}_{}{}_FMAX_results.txt".format(O, self.taxon, Type, Mode, self.author, self.model), 'w')
        FMAXhandle.write('!AUTHOR:      \t {} \n'.format(self.author))
        FMAXhandle.write('!MODEL:       \t {} \n'.format(self.model))
        FMAXhandle.write('!KEYWORDS:    \t {} \n'.format(self.keywords))
        FMAXhandle.write('!SPECIES:     \t {} \n'.format(self.taxon))
        FMAXhandle.write('!ONTOLOGY:    \t {} \n'.format(O))
        FMAXhandle.write('!TYPE:        \t {} \n'.format(Type))
        FMAXhandle.write('!MODE:        \t {} \n'.format(Mode))
        FMAXhandle.write('<OPTIMAL:\t {:.6f} \n'.format(self.FMAX))
        FMAXhandle.write('<THRESHOLD:\t {:.2f} \n'.format(self.FMAXThreshold))
        FMAXhandle.write('<COVERAGE:\t {:.2f} \n'.format(coverage))

        
        FMAXhandle.write('{}\t {}\t {}\t {}\n'.format('Threshold', 'PR', 'RC', 'F')) 
        index = 0
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            try:
                FMAXhandle.write('>{:.2f}\t {:.6f}\t {:.6f}\t {:.6f}\n'.format(threshold, self.PR[index], self.RC[index], self.F[index]))
            except IndexError:
                pass
            index += 1
        
        FMAXhandle.close()
        
        WFMAXhandle  = open(p + "/{}_{}_{}_{}_{}{}_WFMAX_results.txt".format(O, self.taxon, Type, Mode, self.author, self.model), 'w')
        WFMAXhandle.write('!AUTHOR:      \t {} \n'.format(self.author))
        WFMAXhandle.write('!MODEL:       \t {} \n'.format(self.model))
        WFMAXhandle.write('!KEYWORDS:    \t {} \n'.format(self.keywords))
        WFMAXhandle.write('!SPECIES:     \t {} \n'.format(self.taxon))
        WFMAXhandle.write('!ONTOLOGY:    \t {} \n'.format(O))
        WFMAXhandle.write('!TYPE:        \t {} \n'.format(Type))
        WFMAXhandle.write('!MODE:        \t {} \n'.format(Mode))
        WFMAXhandle.write('<OPTIMAL:\t {:.6f} \n'.format(self.WFMAX))
        WFMAXhandle.write('<THRESHOLD:\t {:.2f} \n'.format(self.WFMAXThreshold))
        WFMAXhandle.write('<COVERAGE:\t {:.2f} \n'.format(coverage))        
        

        WFMAXhandle.write('{}\t {}\t {}\t {}\n'.format('Threshold', 'WPR', 'WRC', 'WF'))
        index = 0
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            try:
                WFMAXhandle.write('>{:.2f}\t {:.6f}\t {:.6f}\t {:.6f}\n'.format(threshold  , self.WPR[index], self.WRC[index], self.WF[index]))
            except IndexError:
                pass
            index += 1
        
        WFMAXhandle.close()
        
        SMINhandle  = open(p + "/{}_{}_{}_{}_{}{}_SMIN_results.txt".format(O, self.taxon, Type, Mode, self.author, self.model), 'w')
        SMINhandle.write('!AUTHOR:      \t {} \n'.format(self.author))
        SMINhandle.write('!MODEL:       \t {} \n'.format(self.model))
        SMINhandle.write('!KEYWORDS:    \t {} \n'.format(self.keywords))
        SMINhandle.write('!SPECIES:     \t {} \n'.format(self.taxon))
        SMINhandle.write('!ONTOLOGY:    \t {} \n'.format(O))
        SMINhandle.write('!TYPE:        \t {} \n'.format(Type))
        SMINhandle.write('!MODE:        \t {} \n'.format(Mode))
        SMINhandle.write('<OPTIMAL:\t {:.6f} \n'.format(self.SMIN))
        SMINhandle.write('<THRESHOLD:\t {:.2f} \n'.format(self.SMINThreshold))
        SMINhandle.write('<COVERAGE:\t {:.2f} \n'.format(coverage))          
        
        SMINhandle.write('{}\t {}\t {}\t {}\n'.format('Threshold', 'RU', 'MI', 'S'))
        index = 0
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            try:
                SMINhandle.write('>{:.2f}\t {:.6f}\t {:.6f}\t {:.6f}\n'.format(threshold  , self.RU[index], self.MI[index], self.S[index]))
            except IndexError:
                pass
            index += 1
        
        SMINhandle.close()
        
        NSMINhandle  = open(p + "/{}_{}_{}_{}_{}{}_NSMIN_results.txt".format(O, self.taxon, Type, Mode, self.author, self.model), 'w')
        NSMINhandle.write('!AUTHOR:      \t {} \n'.format(self.author))
        NSMINhandle.write('!MODEL:       \t {} \n'.format(self.model))
        NSMINhandle.write('!KEYWORDS:    \t {} \n'.format(self.keywords))
        NSMINhandle.write('!SPECIES:     \t {} \n'.format(self.taxon))
        NSMINhandle.write('!ONTOLOGY:    \t {} \n'.format(O))
        NSMINhandle.write('!TYPE:        \t {} \n'.format(Type))
        NSMINhandle.write('!MODE:        \t {} \n'.format(Mode))
        NSMINhandle.write('<OPTIMAL:\t {:.6f} \n'.format(self.NSMIN))
        NSMINhandle.write('<THRESHOLD:\t {:.2f} \n'.format(self.NSMINThreshold))
        NSMINhandle.write('<COVERAGE:\t {:.2f} \n'.format(coverage))          
        
        NSMINhandle.write('{}\t {}\t {}\t {}\n'.format('Threshold', 'NRU', 'NMI', 'NS' ))
        index = 0
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            try:
                NSMINhandle.write('>{:.2f}\t {:.6f}\t {:.6f}\t {:.6f}\n'.format(threshold  , self.NRU[index], self.NMI[index], self.NS[index]))
            except IndexError:
                pass
            index += 1
        
        NSMINhandle.close()        
        
        