import numpy
from assessment_new.Tools import vprint, vwrite, clear


def AUC(Info):
    '''
    
    '''

    # Calculate base values
    for term in Info.termPrediction:    
        # Set local copies
        truth      = Info.termTruth[term]
        terms      = Info.termlist[term]
        prediction = Info.termPrediction[term]
        
        # This is the union of terms and truths 
        totalTerms = terms | truth
        # true is all the terms in truth 
        true = truth
        vprint("Positive Set: {}".format(true), 15)
        # Set difference
        # false is the terms that are not true
        false = terms - truth
        vprint("Negative Set: {}".format(false), 15)
        #####################################
        # Set paths
        path = Info.ResultPath + "/Term/{}_{}".format((Info.ontology.lower()), Info.Type)
        clear("{}_{}.txt".format(path, term))      
        message = "Terms: {} \n Truth: {} \n Prediction: {}\n".format(
        terms, truth, prediction)
        vwrite(message, "{}_{}.txt".format(path, term), 1)
        ######################################
        
        # Get counts of terms, constant for all thresholds
        # True (TP + FN)
        TRUE = len(true)
        # False (FP + TN)
        FALSE = len(false) 
################## Calculate FPR & TPR @ all thresholds #######################
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            threshold = numpy.around(threshold, decimals = 2)
            # Intialize Values
            FP  = 0
            TN  = 0
            TP  = 0
            FN  = 0 
            # Get FP, TP values   
            for term in prediction:
                # If Predicted but False
                if (prediction[term] >= threshold 
                    and term in false):
                    FP += 1
                # If Predicted and True
                if (prediction[term] >= threshold 
                    and term in true):
                    TP += 1       
            # Set all values            
            # False Positive
            FP  = FP
            # True Negative
            TN  = FALSE - FP
            # True Positve
            TP  = TP
            # False Negative
            FN  = TRUE - TP
            # Positive Prediction
            POS = TP + FP
            # Negative Prediction (Not Predicted)
            NEG = FN + TN
            # True positive rate (TP / (TP + FN))
            TPR = TP / TRUE
            # False positive rate (FP / (FP + TN))
            FPR = FP / FALSE
            #SPECIAL CASE T = 0 
            if(threshold == 0.00):                
                TP = TRUE                
                POS = Info.OBOCount
            
            # Data_unweighted is for FMAX (Not using IC values)
            Info.data_terms[term][threshold] = {
             'FP':FP, 'TN':TN, 'TP':TP, 'FN':FN, 
             'POS':POS, 'NEG':NEG, 'TRUE':TRUE,'FALSE':FALSE,
             'TPR':TPR, 'FPR':FPR}
            ################
            # Split for readablity 
            part1 = "\t FP: {},\t TN: {},\t TP: {},\t FN: {}".format(
            FP, TN, TP, FN)
            part2 = "\t POS: {},\t NEG: {},\t TRUE: {},\t FALSE: {}".format(
            POS, NEG, TRUE, FALSE)
            part3 = "\t TPR: {},\t FPR: {}".format(TPR, FPR)
            name = "{} @ {:.2f} :".format(term, threshold)
            vwrite("{}: {}, {}, {}\n".format(name, part1, part2, part3), "{}_unweighted.txt".format(path), 1)

    # AUC Caluclations
    # Take data in Info.data_terms to construct areas and return largest area and its threshold
