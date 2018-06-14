import numpy
from assessment.Tools import vprint
#FMAX (send weighted value to use for WFMAX)
def FMAX(Info):
    Fmax          = 0
    FmaxThreshold = -1
    
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
        threshold = numpy.around(threshold, decimals = 2)
        vprint("Threshold is {}".format(threshold),5)
        pr = PR(Info, threshold)
        rc = RC(Info, threshold)
        Fval = F(pr,rc)
        if Fval > Fmax:
            Fmax = Fval
            FmaxThreshold = threshold
    return [Fmax, FmaxThreshold]
    
def PR(Info, threshold):
    # Sum of all proteins 
    total = 0
    for protein in Info.prediction:
        try:
            # TP
            numerator = Info.data_unweighted[protein][threshold]['TP']
            # TP + FP
            demoninator = Info.data_unweighted[protein][threshold]['POS']
            try:
                total += numerator / demoninator
            except ZeroDivisionError:
                vprint("Protein {} had a div by 0".format(protein),5)
        except KeyError:
            vprint("Protein: {}".format(protein),5)
        
    precision = total / Info.ProteinInPrediction[threshold] 
    return precision

def RC(Info, threshold):
    # Sum of all proteins 
    total = 0
    for protein in Info.prediction:
        try:
            # TP
            numerator = Info.data_unweighted[protein][threshold]['TP']
            # TP + FN
            demoninator = Info.data_unweighted[protein][threshold]['TRUE']  
            try:
                total += numerator / demoninator
            except ZeroDivisionError:
                vprint("Protein {} had a div by 0".format(protein),5)
        except KeyError:
            vprint("Protein: {}".format(protein),5)
            
    if Info.mode == "full":
        recall = total / Info.ProteinInBenchmark
    else: # "partial"
        recall = total / Info.ProteinInPrediction[0.00]     
    return recall

def F(pr, rc):
    try:
        return (2*pr*rc)/(pr + rc)
    except ZeroDivisionError:
        return 0