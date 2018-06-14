import math
import numpy
from assessment.Tools import vprint

#SMIN (send weighted value to use for WFMAX)
def SMIN(Info):
    Smin          = float("inf")
    SminThreshold = -1
    
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
        threshold = numpy.around(threshold, decimals = 2)
        vprint("Threshold is {}".format(threshold),5)
        ru = RU(Info, threshold)
        mi = MI(Info, threshold)
        Sval = S(ru, mi)
        if Sval < Smin:
            Smin = Sval
            SminThreshold = threshold
    return [Smin, SminThreshold]
    
def RU(Info, threshold):
    # Sum of all proteins 
    total = 0
    for protein in Info.prediction:
        try:
            total += Info.data[protein][threshold]['FN']
        except KeyError:
            vprint("Protein: {}".format(protein),5)
    if Info.mode == "full":
        remain = total / Info.ProteinInBenchmark
    else: # "partial"
        remain = total / Info.ProteinInPrediction[0.00]    
    return remain

def MI(Info, threshold):
    # Sum of all proteins 
    total = 0
    for protein in Info.prediction:
        try:
            total += Info.data[protein][threshold]['FP']
        except KeyError:
            vprint("Protein: {}".format(protein),5)
    if Info.mode == "full":
        misinfo = total / Info.ProteinInBenchmark
    else: # "partial"
        misinfo = total / Info.ProteinInPrediction[0.00]     
    return misinfo

def S(ru, mi):
    return math.sqrt(ru**2 + mi**2)