import math
import numpy
from assessment_new.Tools import vprint, vwrite, clear


def NSMIN(Info):
    '''
    Calculate the normalized minimum semantic distance for a CAFA prediction
    
    Input:
    Info : Object
    
    Output:
    [0]  : List    [Float, Float]
    '''
    
    # Intialize paths
    path = Info.ResultPath + "/NSMIN/{}/{}/{}".format((Info.ontology.lower()), Info.Type, Info.mode)
    overview = "{}/NSMIN_Overview.txt".format(path)
    clear(overview)
    # Intilize S-val, Threshold
    Smin          = float("inf")
    SminThreshold = -1
    # For all thresholds
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
        threshold = numpy.around(threshold, decimals = 2)
        # Set path for this threshold
        data = "{}/{}/NSMIN_Data.txt".format(path, threshold)
        # Delete old threshold file
        clear(data)
        # Store for inner methods
        Info.local_path = data       
        vprint("Threshold is {}".format(threshold), 15)
        # RU for this prediction @ threshold
        ru = RU(Info, threshold)
        vprint("RU is {}".format(ru), 15)
        # MI for this prediction @ threshold
        mi = MI(Info, threshold)
        vprint("MI is {}".format(mi), 15)
        # S-val for this prediction @ threshold
        Sval = S(ru, mi)
        vprint("The NS-val at {:.2f} is {}".format(threshold, Sval), 15)
        # Write to overview
        vwrite("Threshold: {:.2f},\t RU: {:.4f},\t MI: {:.4f},\t NS: {:.4f}\n".format(threshold, ru, mi, Sval), overview, 1)
        # Check if S-val is less than current S-min
        if Sval < Smin:
            Smin = Sval
            SminThreshold = threshold
    # Clear local_path when done
    Info.local_path = ""
    # Return the S-min and its threshold
    return [Smin, SminThreshold]
    
    
def RU(Info, threshold):
    '''
    Calculates the RU (Remaining Uncertainity) for a given prediction at a threshold
    
    Input:
    Info      : Object     
    threshold : Float       {0.00 -> 1.00}
    
    Output:
    [0]       : Float       Remaining Uncertainity
    '''
    
    # Sum of all proteins 
    total = 0
    for protein in Info.prediction:
        try:
            # Get the FN (False Negatives)
            FN  = Info.data[protein][threshold]['FN']
            # Get the TOT (Predicted or Truth)
            TOT = Info.data[protein][threshold]['TOT']
            vwrite("Protein: {},\t FN: {:.2f},\t TOT: {:.2f}\n".format(protein, FN, TOT), Info.local_path, 1)
            try:
                total += FN / TOT
            except ZeroDivisionError:
                vprint("Protein {} had 0 TOT".format(protein), 4)
        # Bad Error
        except KeyError:
            vprint("Protein: {} has no TOT".format(protein), 1)
    # If in full mode, divide by N (# proteins in benchmark)
    if Info.mode == "full":
        div = Info.ProteinInBenchmark
        remain = total / div
    # If in partial mode, divide by m(0) (# proteins in prediction)
    else: # "partial"
        div = Info.ProteinInPrediction[0.00]
        remain = total / div
    vwrite("Total: {:.2f},\t Ne: {},\t RU: {:.2f}\n".format(total, div, remain) ,Info.local_path, 1)
    # Return the calculated Remaining Uncertainity
    return remain

def MI(Info, threshold):
    '''
    Calculates the MI (Misinformation) for a given prediction at a threshold
    
    Input:
    Info      : Object     
    threshold : Float       {0.00 -> 1.00}
    
    Output:
    [0]       : Float       Misinformation
    '''
    
    # Sum of all proteins 
    total = 0
    for protein in Info.prediction:
        try:
            # Get the FP (False Positive)
            FP = Info.data[protein][threshold]['FP']
            # Get the TOT (Predicted or Truth)
            TOT = Info.data[protein][threshold]['TOT']
            vwrite("Protein: {},\t FP: {:.2f},\t TOT: {:.2f}\n".format(protein, FP, TOT), Info.local_path, 1)
            try:
                total += (FP / TOT)
            # TOT was 0
            except ZeroDivisionError:
                vprint("Protein {} had 0 TOT".format(protein), 4)
        # Bad Error        
        except KeyError:
            vprint("Protein: {}".format(protein),5)
    # If in full mode, divide by N (# proteins in benchmark)
    if Info.mode == "full":
        div = Info.ProteinInBenchmark
        misinfo = total / div
    # If in partial mode, divide by m(0) (# proteins in prediction)
    else: # "partial"
        div = Info.ProteinInPrediction[0.00]  
        misinfo = total / div
    vwrite("Total: {:.2f},\t Ne: {},\t MI: {:.2f}\n".format(total, div, misinfo), Info.local_path, 1)
    # Return the calculated Misinformation    
    return misinfo
    
    
def S(ru, mi):
    '''
    The semanitic distance formula
    
    Input:
    ru  : Float     RemainingUncertainity
    mi  : Float     Misinformation
    
    Output:
    [0] : Float     S-value
    '''
    
    return math.sqrt(ru**2 + mi**2)