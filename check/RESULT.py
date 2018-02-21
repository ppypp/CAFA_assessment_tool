# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 13:48:31 2018

@author: mcgerten
"""

class result:
    ''' Stores results in a common format '''   
    
    
    def __init(self):
        ''' State all variables needed '''
    
        #Fmax
        self.FMAX           = 0.0
        self.PR             = []
        self.RC             = []
        self.FMAXThreshold  = 0.0
        
        #
        self.WFMAX          = 0.0
        self.WPR            = []
        self.WRC            = []
        self.WFMAXThreshold = 0.0
        
        self.SMIN           = 0.0
        self.RU             = []
        self.MI             = []
        self.SMINThreshold  = 0.0
        
        self.NSMIN          = 0.0
        self.NRU            = []
        self.NMI            = []
        self.NSMINThreshold = 0.0
        
        
    def update(self, type, value, subval1, subval2, threshold):
        '''
        Stores the outpu of each metric in a single object
        
        Input:
        type      : String
        value     : Float
        subval1   : List
        subval2   : List
        threshold : Float
        '''
        if type == "FMAX":
            self.FMAX           = value
            self.PR             = subval1
            self.RC             = subval2
            self.FMAXThreshold  = threshold
        
        elif type == "WFMAX":
            self.WFMAX          = value
            self.WPR            = subval1
            self.WRC            = subval2
            self.WFMAXThreshold = threshold
        elif type == "SMIN":
            self.SMIN           = value
            self.RU             = subval1
            self.MI             = subval2
            self.SMINThreshold  = threshold           
        elif type == "NSMIN":
            self.NSMIN          = value
            self.NRU            = subval1
            self.NMI            = subval2
            self.NSMINThreshold = threshold     
        else: 
            # Do Nothing
            return
            
    def printOut(self):
        '''
        
        '''
        print('FMAX: %s\n' % self.FMAX)
        print('threshold giving FMAX: %s\n' % self.FMAXThreshold)
        print('WFMAX: %s\n' % self.WFMAX)
        print('threshold giving WFMAX: %s\n' % self.WFMAXThreshold)
        print('SMIN: %s\n' % self.SMIN)
        print('threshold giving SMIN: %s\n' % self.SMINThreshold)
        print('NSMIN: %s\n' % self.NSMIN)
        print('threshold giving NSMIN: %s\n' % self.NSMINThreshold)
        
        
    def writeOut(self):
        return