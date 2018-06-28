from collections import defaultdict
import numpy
import pickle as cp
import sys
import os
import helper
from Ontology.IO import OboIO
import time 


# Define constants
legal_species = [
 "all", "eukarya", "prokarya",
 'HELPY', 'ECOLI', 'RAT',
 'DANRE', 'SULSO', 'DROME',
 'PSEPK', 'STRPN', 'PSEAE',
 'BACSU', 'MYCGE', 'HUMAN',
 'METJA', 'DICDI', 'YEAST',
 'SCHPO', 'ARATH', 'XENLA',
 'MOUSE', 'PSESM', 'SALTY',
 'CANAX', 'SALCH']
root_terms      = ['GO:0008150', 'GO:0005575', 'GO:0003674']
legal_ontology  = ["BPO", "MFO", "CCO", "HPO"] 
legal_types     = ["type1", "type2", "typex"]
legal_subtypes  = ["easy", "hard"]
legal_modes     = ["full", "partial"]


class Info:
    def __init__(self, ontology, Type):
        '''
        Intialize variables 
        '''
        
        # KEY: Protein              Value: GoTerms
        self.terms                  = defaultdict(set)
        # KEY: Protein              Value: GoTerms
        self.truth                  = defaultdict(set)
        # KEY: Protein              Value: GoTerms
        # KEY: GoTerm               Value: Confidence
        self.prediction             = defaultdict(defaultdict)
        # KEY: Protein              Value: MAX Confidence
        self.maxConfidence          = defaultdict()
        # KEY: GoTerm               Value: IC Value
        self.IC                     = defaultdict()
        # KEY: Protein              Value: Thresholds  
        # KEY: Threshold            Value: {FP, TN, TP, FN, POS, NEG, TRUE, FALSE, TOT}
        self.data                   = defaultdict(defaultdict)
        # KEY: Protein              Value: Thresholds   
        # KEY: Threshold            Value: {FP, TN, TP, FN, POS, NEG, TRUE, FALSE}
        self.data_unweighted        = defaultdict(defaultdict)
        # KEY: threshold            Value: Count of Proteins 
        self.ProteinInPrediction    = defaultdict()
        self.ProteinInBenchmark     = 0
        # FOR AUTHOR / FILE INFO
        self.author                 = ""
        self.model                  = ""
        self.keywords               = ""
        self.taxon                  = ""
    
        # FOR Choosing metric / data
        self.ontology               = ontology
        self.Type                   = Type
        self.mode                   = ""
        
        self.Coverage               = -1
        
        # FOR Benchmark inporting
        self.Benchmark              = None
        self.OBOCount               = 0
        self.OBOICTotal             = 0
        # KEY: GoTerm                  Value: GoTerm
        self.TermAncestors          = defaultdict()
        # KEY: Protein                 Value: GoTerm
        self.ProteinTrueTerms       = defaultdict()
        self.TrueTermCount          = 0
        
        # For filepaths
        self.ResultPath             = ""
        self.PredictionPath         = ""
        self.OBOPath                = ""
        self.Path                   = ""
        # For writing to file, saved here so i dont have to pass the path around
        self.local_path             = ""
        
        # AUC 
        # KEY: GoTerm              Value: Proteins
        self.termList  = defaultdict(set)
        # KEY: GoTerm              Value: Proteins
        self.termTruth = defaultdict(set)
        # KEY: GoTerm              Value: Proteins
        # KEY: Protein             Value: Confidence
        self.termPrediction = defaultdict(defaultdict)
        # KEY: GoTerm              Value: Thresholds   
        # KEY: Threshold           Value: {FP, TN, TP, FN, POS, NEG, TRUE, FALSE, TPR, FPR}
        self.data_terms        = defaultdict(defaultdict)
        
        
    def setMode(self, mode):
        '''
        Set the mode in which the metrics are run
        Don't need to reset Benchamrk to change mode
        '''
        
        self.mode = mode
        
        
    def setTeamInfo(self, author, model, keywords, taxon):
        '''
        Set the team info for print / write out in results
        '''
        
        vprint("Saving Team Info", 9)
        self.author   = author
        vprint("Author: {}".format(self.author), 9)
        self.model    = model
        vprint("Model: {}".format(self.model), 9)
        self.keywords = keywords
        vprint("Keywords: {}".format(self.keywords), 9)
        self.taxon    = taxon
        vprint("Taxon: {}".format(self.taxon), 9)
                
                    
    def setBenchmark(self, OBOCounts, benchmark_directory):
        '''
        Read in the benchamrk file for choosen Taxon/Ontology/Type
        
        Input:
        OBOcounts : Dictionary   KEY: String VALUE: Integer
                                 Ontology    Count of OBO
        benchmark_directory : String 
        '''
        
        # Convert to format files are saved in
        taxon = helper.taxon_name_converter(self.taxon)
        # Read in benchmark
        Benchmark = read_benchmark(self.ontology, taxon, self.Type, benchmark_directory, self.OBOPath)
        self.Benchmark        = Benchmark
        vprint("Benchmark Done", 10)
        # Store OBO counts
        self.OBOCount         = OBOCounts[self.ontology]
        vprint("OBO Count for {} : {}".format(self.ontology, self.OBOCount), 10)
        # Store data from benchmark directly in Info
        self.TermAncestors    = Benchmark.ancestors
        self.ProteinTrueTerms = Benchmark.ProteinTrueTerms
        #Set path
        path  = self.ResultPath + "/{}_{}_Ancestors.txt".format(self.ontology, self.Type)
        clear(path)
        for term in self.TermAncestors:        
            vwrite("{} : {}\n".format(term, self.TermAncestors[term]), path, 1)
        #Set path
        path  = self.ResultPath + "/{}_{}_TrueTerms.txt".format(self.ontology, self.Type)
        clear(path)
        for protein in self.ProteinTrueTerms:        
            vwrite("{} : {}\n".format(protein, self.ProteinTrueTerms[protein]), path, 1)
            
        
    def setIC(self, ic_path):
        '''
        Read in the IC file for choosen Ontology
        MUST RUN IC-Tool of choice once before any prediction evaulation
        ############# IC is in NATS (Log e) currently ############
                      Should be in bits (Log 2)
        '''
        
        # Load dictionary using pickle
        IC_old = cp.load(open("{}ia_{}.map".format(ic_path, self.ontology), "rb"))
        IC = {}
        total = 0
        # Convert legacy IC format to format used in the program
        for term in IC_old:
            vprint("Term : {}".format(term), 10 )
            vprint("IC   : {}".format(IC_old[term][1]), 10)
            IC[term] = IC_old[term][1]
            total += IC[term]
        # Store convert IC values
        self.IC = IC   
        self.OBOICTotal = total
        vprint(self.IC, 15)
        #Set path
        path  = self.ResultPath + "/IC.txt"
        clear(path)
        for term in IC:        
            vwrite("{} : {:.2f}\n".format(term, IC[term]), path, 1)
        
        
    def setPaths(self, results_path, prediction_path, obo_path, path):  
        '''
        Set paths for necessary files
        '''
        
        self.ResultPath = results_path
        vprint("Results Path: {}".format(self.ResultPath), 9)
        # Make directory for results / data
        helper.mkdir_results(self.ResultPath)
        vprint("Result directory structure complete", 9)
        self.PredictionPath = prediction_path
        vprint("Prediction Path: {}".format(self.PredictionPath), 9)
        self.OBOPath = obo_path
        vprint("OBO Path: {}".format(self.OBOPath), 9)
        self.Path = path
        vprint("Path: {}".format(self.Path), 9)
        
        
    def setPrediction(self):
        '''
        Read in the prediction file and convert to usable format
        '''
        
        # Intialize variables
        self.ProteinInBenchmark = 0
        PredictionProteinTemp = defaultdict(list)  
        # If the file exists
        if os.path.isfile(self.Path):
            # Read in prediction file   
            handle = open(self.Path,'r')
            for inrec in handle:
                fields     = [i.strip() for i in inrec.split()]
                protein    = fields[0]
                vprint("Protein: {}".format(protein), 15)
                term       = fields[1]
                vprint("Term: {}".format(term), 15)
                confidence = float(fields[2])
                vprint("Confidence: {}".format(confidence), 15)
                # Store prediction in temporary dictionary
                PredictionProteinTemp[protein].append(
                        {'term': term, 'confidence': confidence})
            # Close to be filesafe
            handle.close()
            # Propagate prediction
            for protein in PredictionProteinTemp:
                vprint("Protein: {}".format(protein), 15)
                # Take only Protein that has terms in the benchmark
                if self.ProteinTrueTerms[protein]:
                    self.ProteinInBenchmark += 1
                    vprint(PredictionProteinTemp[protein], 15)
                    # For all pairs in PredictionProteinTemp[protein], pull the data
                    for value in PredictionProteinTemp[protein]:
                        term       = value['term']
                        confidence = value['confidence']
                        # Get term ancestors                         
                        try:
                            ancterms = self.TermAncestors[term].difference(root_terms)
                            vprint(ancterms,15)
                        except KeyError:
                            # Add unknown term to the obsolete set
                            self.obsolete.add(term)
                            continue
                        # For each term
                        if term in self.prediction[protein]:
                            # Term already exists, update confidence
                            self.update_confidence(protein, term, confidence)
                        else:
                            # Add term to self.predicted_bench
                            if term not in root_terms:
                                self.prediction[protein][term] = confidence
                                # Propagate for ancestors
                                for ancterm in ancterms:
                                    if ancterm in self.prediction[protein]:
                                        # Term already exists, update confidence
                                        self.update_confidence(protein, term, confidence)
                                    else:
                                        # Add term to self.predicted_bench
                                        self.prediction[protein][ancterm] = confidence            
            if self.ProteinInBenchmark == 0:
                vprint("No protein in this predicted set became a benchmark\n", 1)
            else:
                vprint("Protein in Benchmark: {}".format(self.ProteinInBenchmark), 10)
        else:
            vprint('No prediction made in this ontology.\n', 1)
           
            
    def coverage(self):
        ''' 
        Determine the coverage ( m(0) / N ). 
        
        Output:
        [0] : Float       coverage of the prediction
        '''
        
        try:
            cover = float(self.ProteinInPrediction[0.00])/self.ProteinInBenchmark
            return cover
        except ZeroDivisionError:
            return None
        
        
    def update_confidence(self, protein, term, confidence):
        '''
        Update Confidence for given protein and propagate.
        
        This function compares the confidence value given to the confidence in self.predicted
        If the given is larger, then it overwrites the confidence in self.predicted
        
        Input:
        protein    : String
        term       : String
        confidence : Float
        '''
        
        if confidence > self.prediction[protein][term]:
            # Update the confidence
            self.prediction[protein][term] = confidence
            # Propagate changes if necessary
            for ancterm in self.TermAncestors[term].difference(root_terms):
                if confidence > self.prediction[protein][ancterm]:
                    # Update the confidence
                    self.prediction[protein][ancterm] = confidence
                            
        
    def getObsolete(self):
        ''' 
        Get obsolete terms used by the prediction team. 
        
        Output:
        [0] : Set      GOTerms that are obsolute
        '''
        
        return(self.obsolete)   
    
    
    def setInfo(self):
        '''
        Turn the inputed prediction and benchmark data 
        into a format usable by the metrics
        This approach seemed to make more sense as it calculates everything 
        once rather than each time a metric is called
        '''
        
        # Find Max confidence for each protein
        # We need this in order to find m(T) 
        for protein in self.prediction:
            self.maxConfidence[protein] = 0
            for term in self.prediction[protein]:
                if self.prediction[protein][term] > self.maxConfidence[protein]:
                    self.maxConfidence[protein] = self.prediction[protein][term]
        vprint("Max Confidence : {} ".format(self.maxConfidence), 15)    
        
        # Generate protein in prediction count for each threshold ( m(T) )
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            threshold = numpy.around(threshold, decimals = 2)
            total = 0
            # Check if protein is above/equal to the is threshold
            for protein in self.prediction:
                if self.maxConfidence[protein] >= threshold:
                    total += 1
            # Store the m(T) values
            self.ProteinInPrediction[threshold] = total
            vprint("ProteinInPrediction at {} : {}".format(threshold, total), 10)
        
        # Generate Terms and truths
        for protein in self.prediction:
            terms  = set()
            truths = self.ProteinTrueTerms[protein]
            # Grab the terms from the dictionary
            for term in self.prediction[protein]:
                terms.add(term)
            self.terms[protein] = terms
            self.truth[protein] = truths
            vprint(protein, 15)
            vprint("Terms: {}".format(terms), 15)
            vprint("Truth: {}".format(truths), 15)
        
        print("Info has been set")
        # Generate other need values
        self.generateIntermediateData()
        # Run coverage 
        self.Coverage = self.coverage()
        # Run AUC Benchmark
        self.convertbenchmarkforAUC()
        # Ready to run metrics now
        return True 
    
    
    ################### Actual work for metrics ##################   
    def generateIntermediateData(self):
        '''
        Generate and store the stat values in order to calculate metrics
        '''
        
        # Set paths
        path = self.ResultPath + "/{}_{}".format((self.ontology.lower()), self.Type)
        clear("{}.txt".format(path))
        clear("{}_unweighted.txt".format(path))
        # For all proteins in prediction
        for protein in self.prediction:
            # Set local variables
            terms       = self.terms[protein]
            truth       = self.truth[protein]
            prediction  = self.prediction[protein]
            ic          = self.IC
            # This is the union of terms and truths 
            totalTerms = terms | truth
            '''# Get the subset of terms that have IC values
            terms_with_IC = terms & ic.keys()
            # Create a list of terms that dont have IC values 
            badTerms = terms - terms_with_IC
            vprint("Bad terms Set: {}".format(badTerms), 15)
            # Bad terms are those without IC
            # Remove bad terms from terms
            terms = terms - badTerms
            '''
            ########################################################
            '''# ????
            # Remove IC-less terms from prediction   
            # Get Prediction Terms that have IC
            prediction_with_IC = prediction.keys() & ic.keys()
            # Then find those that dont have IC
            badTerms = prediction.keys() - prediction_with_IC
            for term in badTerms:
                # Remove term that doesn't have IC
                try:                
                    prediction.pop(term)
                # The term was not in the prediction
                except KeyError:
                    pass
            # NOW the prediction dictionary won't have non-IC terms
            '''# IS this desired, not sure
            ########################################################
            vprint("Terms Set: {}".format(terms), 15)
            # true is all the terms in truth (ProteinTrueTerms) 
            true = truth
            vprint("Positive Set: {}".format(true), 15)
            # Set difference
            # false is the terms that are not true
            false = terms - truth
            vprint("Negative Set: {}".format(false), 15)
            #####################################
            # Set paths
            path = self.ResultPath + "/Protein/{}_{}".format((self.ontology.lower()), self.Type)
            clear("{}_{}.txt".format(path, protein))      
            message = "Terms: {} \n Truth: {} \n Prediction: {}\n".format(
            terms, truth, prediction)
            vwrite(message, "{}_{}.txt".format(path, protein), 1)
            ######################################
            
            # Get counts of terms, constant for all thresholds
            # True (TP + FN)
            TRUE = len(true)
            # False (FP + TN)
            FALSE = len(false) 
            # For each threshold, find values (For FMAX)
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
                
                #SPECIAL CASE T = 0 
                if(threshold == 0.00):                
                    TP = TRUE                
                    POS = self.OBOCount
                
                # Data_unweighted is for FMAX (Not using IC values)
                self.data_unweighted[protein][threshold] = {
                 'FP':FP, 'TN':TN, 'TP':TP, 'FN':FN, 
                 'POS':POS, 'NEG':NEG, 'TRUE':TRUE,'FALSE':FALSE}
                ################
                # Split for readablity 
                part1 = "\t FP: {},\t TN: {},\t TP: {},\t FN: {}".format(
                FP, TN, TP, FN)
                part2 = "\t POS: {},\t NEG: {},\t TRUE: {},\t FALSE: {}".format(
                POS, NEG, TRUE, FALSE)
                name = "{} @ {:.2f} :".format(protein, threshold)
                vwrite("{}: {}, {}\n".format(name, part1, part2), "{}_unweighted.txt".format(path), 1)
                ################
            #COUNT FOR WEIGHTED METRICS
            # Intialize Values
            weightedTrue = 0
            weightedFalse = 0
            # Get weightedTrue, weightedFalse values
            for term in totalTerms:
                # If terms is true, add IC to weightedTrue
                if term in true:
                    try:
                        weightedTrue += ic[term]
                        vprint("{} added to TRUE".format(ic[term]), 7)
                    # Shouldn't happen as we removed terms that are missing IC    
                    except KeyError:
                        vprint("WP : {} has no IC value".format(term), 5) 
                # If terms is false, add IC to weightedFalse
                if term in false:
                    try:
                        weightedFalse += ic[term]   
                        vprint("{} added to FALSE".format(ic[term]), 7)
                    # Shouldn't happen as we removed terms that are missing IC    
                    except KeyError:
                        vprint("WN : {} has no IC value".format(term), 5) 
            # Store in same format as other values
            # TRUE, FALSE not threshold dependent
            TRUE  = weightedTrue
            FALSE = weightedFalse

            # For all thresholds
            for threshold in numpy.arange(0.00, 1.01, 0.01, float):
                threshold = numpy.around(threshold, decimals = 2)
                # Intialize Values
                FP = 0
                TN = 0
                TP = 0
                FN = 0
                # Get FP, TP values   
                for term in prediction:
                    # If Predicted but False
                    if (prediction[term] >= threshold 
                        and term in false):
                        # Add to FP (False Positive)   
                        try:
                            FP += ic[term]
                            vprint("{} added to FP".format(ic[term]), 7)
                        # Shouldn't happen as we removed terms that are missing IC    
                        except KeyError:
                            vprint("FP : {} has no IC value".format(term), 5) 
                    # If Predicted and True
                    if (prediction[term] >= threshold 
                        and term in true):
                        # Add to TP (True Positive)   
                        try:        
                            TP += ic[term]
                            vprint("{} added to TP".format(ic[term]), 7)
                        # Shouldn't happen as we removed terms that are missing IC    
                        except KeyError:
                            vprint("TP : {} has no IC value".format(term), 5) 
                # Set all values            
                # False Positive
                FP  = FP
                vprint("{} at {:.2f}  has FP : {}".format(protein, threshold, FP), 6)
                # True Negative
                TN  = FALSE - FP
                vprint("{} at {:.2f}  has TN : {}".format(protein, threshold, TN), 6)
                # True Positve
                TP  = TP
                vprint("{} at {:.2f}  has TP : {}".format(protein, threshold, TP), 6)
                # False Negative
                FN  = TRUE - TP
                vprint("{} at {:.2f}  has FN : {}".format(protein, threshold, FN), 6)
                # Positive Prediction
                POS = TP + FP
                vprint("{} at {:.2f}  has POS : {}".format(protein, threshold, FN), 6)
                # Negative Prediction (Not Predicted)
                NEG = FN + TN
                vprint("{} at {:.2f}  has NEG : {}".format(protein, threshold, FN), 6)
                # Used for normalizing S-Min                
                # Predicted or Truth
                # Also could be calculated by (POS + NEG) or (TRUE + FALSE)
                # as both those should have the same IC sum
                TOT = Info.tot(ic, truth, prediction)
                vprint("{} at {:.2f}  has TOT : {}".format(protein, threshold, FN), 6)
                
                #SPECIAL CASE T = 0
                if(threshold == 0.00):                 
                    TP  = TRUE
                    POS = self.OBOICTotal
                    
                '''##############################################
                if protein == "T96060019016":
                    print("TRUE: {}".format(TRUE))
                    print(true)
                    print(terms)
                    print(totalTerms)
                '''##############################################
                    
                # Store values
                self.data[protein][threshold] = {
                  'FP':FP, 'TN':TN, 'TP':TP, 'FN':FN, 
                  'POS':POS, 'NEG':NEG, 'TRUE':TRUE,'FALSE':FALSE, 'TOT':TOT}
                #############################
                # Split for readablity 
                part1 = "\t FP: {:.2f},\t TN: {:.2f},\t TP: {:.2f},\t FN: {:.2f}".format(
                FP, TN, TP, FN)
                part2 = "\t POS: {:.2f},\t NEG: {:.2f},\t TOT: {:.2f},\t TRUE: {:.2f},\t FALSE: {:.2f}".format(
                POS, NEG, TOT, TRUE, FALSE)
                name = "{} @ {:.2f} :".format(protein, threshold)
                vwrite("{}: {}, {}\n".format(name, part1, part2), "{}.txt".format(path), 1)
                #############################
                
        print ("Done with Intermediate value creation")  
    
     ################## Convert Benchmark to Term-centric #############
    def convertbenchmarkforAUC(self):
        '''
        Convert Benchmark to Term-based
        For use with AUC
        
        '''
                
        for protein in self.prediction:
            # Local copies for particular protein
            terms       = self.terms[protein]
            truth       = self.truth[protein]
            prediction  = self.prediction[protein]
            # Construct Lists
            for term in terms:
                self.termList[term].add(protein)
            for term in truth:
                self.termTruth[term].add(protein)
            for term in prediction:
                self.termPrediction[term][protein] = prediction[term]
        #################################
        # Only keep those with 10+ proteins
        #################################
        delList = []
        for term in self.termTruth:
            # If less than 10, delete all references
            if (len(self.termTruth[term]) < 10):
                delList.append(term)
        for term in delList:
            try:
                del(self.termList[term])
            except:
                pass
            try:
                del(self.termTruth[term])
            except:
                pass
            try:            
                del(self.termPrediction[term])
            except:
                pass
        # Now have a benchmark set for terms 
        # where every term has at least 10 proteins in the truth
        
        
        
    def tot(IC, T, P):
        '''
        Sum the IC of all terms
        
        Input:
        info       : Object
        T   : Set [ Truth      ]
        P   : Set [ Prediction ]
        
        Output:
        [0] : Float
        '''
        
        #vprint("Finding TOT", 5) 
        # Sum of IC values of terms meeting criteria   
        total = 0.0
        # Sum the IC values of every element in T or P
        # Add all truths
        for term in T:
            try:
                total += IC[term]
            except KeyError:
                vprint("T : {} has no IC value".format(term), 5)  
        # Add predictions that are not truths
        # In order to not double count
        for term in P:
            if term not in T:
                try:
                    total += IC[term]
                except KeyError:
                    vprint("P !T : {} has no IC value".format(term), 5)  
        # Return the final value
        return total
        
        
    def setTime(self, start_time):
        '''
        Set the start time in order to do time ansylis in the Info object
        '''
        
        self.start_time = start_time
        
        
###################### HELPER FUNCTIONS ###################
def vprint(message, priority):
    '''
    Print function with an adjustable verbosity
    Change number for different verbosity
    '''
    
    if priority < 6:
        print(message)
        
        
def clear(location):
    '''
    Delete old results/ files
    '''
    
    try:
        os.remove(location)
    # If the file didn't exist, ignore
    except OSError:
        pass    
    vprint("Old Results have been cleared", 10)
    
    
def vwrite(message, location, priority):
    '''
    Write function with an adjustable verbosity
    Change number for different verbosity
    '''
    
    if priority < 3:
        # Append to file
        data = open(location, 'a+')
        data.write(message)
        data.close()


def getTime(current_time):
    '''
    Facilate Time broadcasts for optimizing
    '''

    if current_time == 0:
        return time.time()
    else:
        vprint((time.time() - current_time), 2)
            
    
class benchmark:
    '''
    Defines a benchmark for comparing the predictions to
    '''
    
    def __init__(self, ancestor_path, benchmark_path):
        '''
        Initialize the benchmark.
        
        Input: 
        benchmark_path : String    ontology specific file location
        ancestor_path  : String    ontology specific file location
        '''
        
        # Key: protein
        # Value: set of benchmark leaf terms
        self.ancestors = defaultdict(set)
        # Read GO ancestors file generated with go_ontology_ancestors_split_write()
        # File format: 
        # go_term <tab> ancestor_1,ancestor_2,..,ancestor_n
        with open(ancestor_path) as ancestors_input:
            for inline in ancestors_input:
                inrec = inline.strip().split('\t')
                term = inrec[0]
                if len(inrec) == 1:
                    self.ancestors[term] = set({})
                else:
                    term_ancestors = inrec[1]
                    self.ancestors[term] = set(term_ancestors.split(','))
        # Build Truth        
        self.ProteinTrueTermsBase = defaultdict(set)
        with open(benchmark_path) as benchmark_input:
            for inline in benchmark_input:
                protein, term = inline.strip().split('\t')
                self.ProteinTrueTermsBase[protein].add(term)
                
                
    def propagateTrueTerms(self):
        '''
        Progate Benchmark terms.
        '''
        
        self.ProteinTrueTerms = defaultdict(set)
        # Key: protein
        # Value: set of benchmark propagated terms
        
        for protein in self.ProteinTrueTermsBase:
            for term in self.ProteinTrueTermsBase[protein]:
                try:
                    ancestors = self.ancestors[term].difference(root_terms)    
                except KeyError:
                    sys.stderr.write("{} not found \n".format(term)) 
                # Add the term    
                self.ProteinTrueTerms[protein].add(term)
                # Set Union 
                self.ProteinTrueTerms[protein] |= ancestors


def read_benchmark(ontology, species, types, benchmark_directory, obopath):
    '''
    Read Benchmark.
    
    Input:
    ontology            : String        {BPO, CCO, MFO, HPO}
    species             : String        List above of accepted options
    types               : String        {full, partial}
    fullbenchmarkfolder : String        Directory
    obopath             : String        File 
    
    Output:
    [0]                 : Benchmark
    [1]                 : obocountDict 
    '''
    
    # Error Checking
    if ontology not in legal_ontology:
        sys.stderr.write("Namespace not accepted, choose from 'BPO', 'CCO', 'MFO' and 'HPO'\n")
    elif (species not in legal_species) and (species not in legal_subtypes):
        sys.stderr.write('Species not accepted')
    elif types not in legal_types:
        sys.stderr.write('Type not accepted, choose from "type1","type2" and "typex"\n')
    else:
        # Create name for finding benchmark file
        matchname = ontology.lower() + '_' + species + '_' + types + '.txt'
        vprint("Matchname: {}".format(matchname), 9)
    # Ontology-specific calculations
    full_benchmark_path = benchmark_directory + '/groundtruth/' + 'leafonly_{}.txt'.format(ontology)
    vprint("full_benchmark_path: {}".format(full_benchmark_path), 9)
    ancestor_path       = os.path.splitext(obopath)[0] + "_ancestors_{}.txt".format(ontology)
    vprint("ancestor_path: {}".format(ancestor_path), 9)
    benchmarkListPath   = benchmark_directory + '/lists/' + matchname
    vprint("benchmarkListPath: {}".format(benchmarkListPath), 9)
    # If the List file exists and has data
    if os.path.isfile(benchmarkListPath) and os.path.getsize(benchmarkListPath) > 0:
        handle    = open(benchmark_directory + '/lists/' + matchname, 'r')
        proteins  = set()
        for line in handle:
            proteins.add(line.strip())
        handle.close()
        tempfilename = 'temp_{}_{}_{}.txt'.format(ontology, species, types)
        tempfile = open(benchmark_directory + '/' + tempfilename , 'w')
        for line in open(full_benchmark_path,'r'):
            protein = line.split('\t')[0]
            if protein in proteins:
                tempfile.write(line)
        tempfile.close()
        bench = benchmark(ancestor_path, tempfile.name)
        bench.propagateTrueTerms()
        #os.remove(tempfile.name)
    else:
        print('Benchmark set is empty.\n')
        bench = None
    return bench 


############### Read OBO ################  
def go_ontology_split(OBO):
    '''
    Split an GO obo file into three ontologies
    
    Input:
    ontology : Ontology File generate by OBOReader
    
    Output:
    [0]      : Set   MFO Terms
    [1]      : Set   BPO Terms
    [2]      : Set   CCO Terms
    '''
    
    #Intialize sets
    mfo_terms = set({})
    bpo_terms = set({})
    cco_terms = set({})
    # Split terms by their ontology
    for term in OBO.get_ids(): # loop over node IDs and alt_id's
        if OBO.namespace[term]   == "biological_process": # BPO, P
            bpo_terms.add(term)
        elif OBO.namespace[term] == "cellular_component": # CCO, C
            cco_terms.add(term)
        elif OBO.namespace[term] == "molecular_function": # MFO, F
            mfo_terms.add(term)
        else:
            raise(ValueError,"{} has no namespace".format(term))
    # Return the sets
    return (mfo_terms, bpo_terms, cco_terms)
 
   
def readOBO(obo_path):
    '''
    Takes an OBO file and writes 3 files with ancestors corresponding to the namespace
    
    Input:
    obo_path : filepath of the obo file
    
    Output:
    [0]      : List [Length of BPO, Length of CCO, Length of MFO]
    '''
    
    OBO = OboIO.OboReader(open(obo_path)).read()
    # Split the OBO into the three namespaces 
    MFO_terms, BPO_terms, CCO_terms = go_ontology_split(OBO)
    # BPO
    obo_BPO_out = open("{}_ancestors_BPO.txt".format(os.path.splitext(obo_path)[0]), "w")
    for term in BPO_terms:
        ancestors = OBO.get_ancestors(term)
        obo_BPO_out.write("%s\t%s\n" % (term, ",".join(ancestors)))
    obo_BPO_out.close()
    # CCO
    obo_CCO_out = open("{}_ancestors_CCO.txt".format(os.path.splitext(obo_path)[0]), "w")
    for term in CCO_terms:
        ancestors = OBO.get_ancestors(term)
        obo_CCO_out.write("%s\t%s\n" % (term, ",".join(ancestors)))
    obo_CCO_out.close()    
    # MFO
    obo_MFO_out = open("{}_ancestors_MFO.txt".format(os.path.splitext(obo_path)[0]), "w") 
    for term in MFO_terms:
        ancestors = OBO.get_ancestors(term)
        obo_MFO_out.write("%s\t%s\n" % (term, ",".join(ancestors)))    
    obo_MFO_out.close()
    
    OBOcountDict = {'BPO':len(BPO_terms),'CCO':len(CCO_terms),'MFO':len(MFO_terms)}
    return OBOcountDict