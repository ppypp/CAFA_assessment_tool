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
        vprint("PR is {}".format(pr),5)
        rc = RC(Info, threshold)
        vprint("RC is {}".format(rc),5)
        Fval = F(pr,rc)
        vprint("The F-val at {} is {}".format(threshold, Fval),1)
        if Fval >= Fmax:
            Fmax = Fval
            FmaxThreshold = threshold
    return [Fmax, FmaxThreshold]
    
def PR(Info, threshold):
    # Sum of all proteins 
    total = 0
    for protein in Info.prediction:
        # Summation of Sum(TP)/ Sum(POS)
        try:
            # TP
            numerator = Info.data_unweighted[protein][threshold]['TP']
            vprint("TP : {}".format(numerator),15)
            # TP + FP
            demoninator = Info.data_unweighted[protein][threshold]['POS']
            vprint("POS : {}".format(demoninator),15)
            try:
                total += numerator / demoninator
                #vprint("I succesfully added to total",2)
            except ZeroDivisionError:
                vprint("Protein {} had a 0 POS or TP".format(protein),15)
        except KeyError:
            vprint("Protein: {}".format(protein),5)
    
    vprint("ProteinInPrediction: {}".format(Info.ProteinInPrediction[threshold]),5)    
    # Divide by m(T)    
    try:    
        precision = total / Info.ProteinInPrediction[threshold] 
    except ZeroDivisionError:
        vprint("No protein predicted at {}".format(threshold),5)
        precision = None
    return precision

def RC(Info, threshold):
    # Sum of all proteins 
    total = 0
    for protein in Info.prediction:
        try:
            # TP
            numerator = Info.data_unweighted[protein][threshold]['TP']
            vprint("TP : {}".format(numerator),15)
            # TP + FN
            demoninator = Info.data_unweighted[protein][threshold]['TRUE'] 
            vprint("TRUE : {}".format(demoninator),15)
            try:
                total += numerator / demoninator
            except ZeroDivisionError:
                vprint("Protein {} had a 0 TRUE or TP".format(protein),5)
        except KeyError:
            vprint("Protein: {}".format(protein), 5)
            
    if Info.mode == "full":
        vprint("Recall Divide Full: {}".format(Info.ProteinInBenchmark),5)
        recall = total / Info.ProteinInBenchmark
    else: # "partial"
        vprint("Recall Divide Partial: {}".format(Info.ProteinInPrediction[0.00]),5)
        recall = total / Info.ProteinInPrediction[0.00]     
    return recall

def F(pr, rc):
    try:
        return (2*pr*rc)/(pr + rc)
    except ZeroDivisionError:
        return 0