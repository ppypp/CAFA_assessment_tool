import numpy
from assessment_new.Tools import vprint, vwrite, clear


def FMAX(Info):
    '''
    Calculate the maximun harmonic mean (FMAX) for a CAFA prediction
    
    Input:
    Info : Object
    
    Output:
    [0]  : List    [Float, Float]
    '''
    
    # Intialize paths
    path = Info.ResultPath + "/FMAX/{}/{}/{}".format((Info.ontology.lower()), Info.Type, Info.mode)
    overview = "{}/FMAX_Overview.txt".format(path)
    clear(overview)  
    # Intilize F-val, Threshold
    Fmax          = 0
    FmaxThreshold = -1
    # For all thresholds
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
        threshold = numpy.around(threshold, decimals = 2)
        # Set path for this threshold
        data = "{}/{}/FMAX_Data.txt".format(path, threshold)
        # Delete old threshold file
        clear(data)
        # Store for inner methods
        Info.local_path = data       
        vprint("Threshold is {}".format(threshold), 5)
        # PR for this prediction @ threshold
        pr = PR(Info, threshold)
        vprint("PR is {}".format(pr), 5)
        # RC for this prediction @ threshold
        rc = RC(Info, threshold)
        vprint("RC is {}".format(rc), 5)
        # F-val for this prediction @ threshold
        Fval = F(pr,rc)
        vprint("The F-val at {:.2f} is {}".format(threshold, Fval), 5)
        # Write to overview
        vwrite("Threshold: {:.2f},\t PR: {:.4f},\t RC: {:.4f},\t F: {:.4f}\n".format(threshold, pr, rc, Fval), overview, 1)
        # Check if F-val is greater than current F-max
        if Fval >= Fmax:
            Fmax = Fval
            FmaxThreshold = threshold
    # Clear local_path when done
    Info.local_path = ""
    # Return the F-max and its threshold
    return [Fmax, FmaxThreshold]
    
    
def PR(Info, threshold):
    '''
    Calculates the PR (Precision) for a given prediction at a threshold
    
    Input:
    Info      : Object     
    threshold : Float       {0.00 -> 1.00}
    
    Output:
    [0]       : Float       Precision
    '''
    
    # Sum of all proteins 
    total = 0
    for protein in Info.prediction:
        # Summation of Sum(TP)/ Sum(POS)
        try:
            # Get the TP (True Positive)
            TP = Info.data_unweighted[protein][threshold]['TP']
            vprint("TP : {}".format(TP), 15)
            # Get the POS (TP + FP) (Positive)
            POS = Info.data_unweighted[protein][threshold]['POS']
            vprint("POS : {}".format(POS), 15)
            try:
                total += TP / POS
            except ZeroDivisionError:
                vprint("Protein {} had a 0 POS".format(protein), 4)
        # Bad Error
        except KeyError:
            vprint("Protein: {} has no POS".format(protein), 1)
        vwrite('Protein: {}\t TP: {}\t POS: {}\n'.format(protein, TP, POS), Info.local_path, 1)
    # Divide by m(T) (# of proteins with at least one prediction @ threshold)   
    try:    
        precision = total / Info.ProteinInPrediction[threshold] 
    except ZeroDivisionError:
        vprint("No protein predicted at {}".format(threshold), 4)
        precision = 0
    vwrite('Precision: {:.2f}\n'.format(precision), Info.local_path, 1)
    # Return the calculated Precision
    return precision


def RC(Info, threshold):
    '''
    Calculates the RC (Recall) for a given prediction at a threshold
    
    Input:
    Info      : Object     
    threshold : Float       {0.00 -> 1.00}
    
    Output:
    [0]       : Float       Recall
    '''
    # Sum of all proteins 
    total = 0
    for protein in Info.prediction:
        # Summation of Sum(TP)/ Sum(TRUE)
        try:
            # Get the TP (True Positive)
            TP = Info.data_unweighted[protein][threshold]['TP']
            vprint("TP : {}".format(TP),15)
            # Get the TRUE (TP + FN) (Truth)
            TRUE = Info.data_unweighted[protein][threshold]['TRUE'] 
            vprint("TRUE : {}".format(TRUE),15)
            try:
                total += TP / TRUE
            except ZeroDivisionError:
                vprint("Protein {} had a 0 TRUE".format(protein), 15)
        # Bad Error
        except KeyError:
            vprint("Protein: {} has no TRUE".format(protein), 1)
        vwrite('Protein: {}\t TP: {}\t TRUE: {}\n'.format(protein, TP, TRUE), Info.local_path, 1)
    # If in full mode, divide by N (# proteins in benchmark)    
    if Info.mode == "full":
        recall = total / Info.ProteinInBenchmark
    # If in partial mode, divide by m(0) (# proteins in prediction)
    else: # "partial"
        recall = total / Info.ProteinInPrediction[0.00]
    vwrite('Recall: {:.2f}\n'.format(recall), Info.local_path, 1)
    # Return the calculated Recall
    return recall


def F(pr, rc):
    '''
    The harmonic mean formula
    
    Input:
    pr  : Float     Precision
    rc  : Float     Recall
    
    Output:
    [0] : Float     F-value
    '''
    
    try:
        return (2*pr*rc)/(pr + rc)
    except ZeroDivisionError:
        return 0