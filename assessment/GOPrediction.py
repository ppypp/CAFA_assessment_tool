#!/usr/bin/env python3

import re
import os
from collections import defaultdict
from Ontology.IO import OboIO


###################### What does this do #####################################
pr_field         = re.compile("^PR=[0,1]\.[0-9][0-9];$")
rc_field         = re.compile("^RC=[0,1]\.[0-9][0-9]$")
# Fix to add HP ontology 2014-1-9
go_field         = re.compile("^(GO|HP):[0-9]{5,7}$")
# Fix to add EFI targets 2014-1-9
target_field     = re.compile("^(T|EFI)[0-9]{5,20}$")
confidence_field = re.compile("^[0,1]\.[0-9][0-9]$")

# Legal states: the CAFA prediction records fields, and their order. KEYWORDS and ACCURACY are
# optional
legal_states1 = ["author", "model", "keywords",
                 "accuracy", "go_prediction", "end"]
legal_states2 = ["author", "model", "keywords", "go_prediction", "end"]
legal_states3 = ["author", "model", "go_prediction", "end"]
legal_states4 = ["go_prediction"]

legal_keywords = [
    "sequence alignment", "sequence-profile alignment", "profile-profile alignment", "phylogeny",
    "sequence properties",
    "physicochemical properties", "predicted properties", "protein interactions", "gene expression",
    "mass spectrometry",
    "genetic interactions", "protein structure", "literature", "genomic context", "synteny",
    "structure alignment",
    "comparative model", "predicted protein structure", "de novo prediction", "machine learning",
    "genome environment",
    "operon", "ortholog", "paralog", "homolog", "hidden markov model", "clinical data", "genetic data",
    "natural language processing", "other functional information"
]


class GOPrediction:
    """
    A class for reading and storing CAFA GO predictions

    in self.read(pred_path) , pred_path should be a handle
    The read function should read in from a RAW prediction file submitted by a CAFA participanting team (post format check though)
    The split function will perform spliting the prediction by ontology, which used to be a separate function in preprocess.py
    """

    def __init__(self):
        '''
        Initilize the prediction
        '''
        
        # get author (teamID) and model from both header and filename
        self.author   = None
        self.model    = None
        # get keywords and accuracy from header
        self.keywords = []
        # split by ontology
        # get taxon from filename
        self.taxon    = None
        # self.data is a dictionary
        # key: protein ID
        # value: [{'term':go_term_1, 'threshold': threshold_1},...,{'term':go_term_n, 'threshold': threshold_n}]
        self.data     = defaultdict(list)
        
        
    def _author_check(self, inrec):
        '''
        Checking the author field
        
        Input:
        inrec  :
        
        Output:
        [0]    : boolean   Correctness
        [1]    : String    Message
        '''
        
        correct = True
        errmsg  = None
        fields  = [i.strip() for i in inrec.split()]
        if len(fields) != 2:
            correct = False
            errmsg  = "AUTHOR: invalid number of fields. Should be 2"
        elif fields[0] != "AUTHOR":
            correct = False
            errmsg  = "AUTHOR: First field should be AUTHOR"
        else:
            self.author = fields[1]
        return correct, errmsg


    def _model_check(self, inrec):
        '''
        Checking the model field
        
        Input:
        inrec  :
        
        Output:
        [0]    : boolean   Correctness
        [1]    : String    Message
        '''
        
        correct = True
        errmsg  = None
        fields  = [i.strip() for i in inrec.split()]
        if len(fields) != 2:
            correct = False
            errmsg  = "MODEL: invalid number of fields. Should be 2"
        elif fields[0] != "MODEL":
            correct = False
            errmsg  = "MODEL: First field should be MODEL"
        elif len(fields[1]) != 1 or not fields[1].isdigit():
            correct = False
            errmsg  = "MODEL: second field should be single digit."
        else:
            self.model = int(fields[1])
        return correct, errmsg


    def _keywords_check(self, inrec):
        '''
        Checking the keywords field
        
        Input:
        inrec  :
        
        Output:
        [0]    : boolean   Correctness
        [1]    : String    Message
        '''
        correct = True
        errmsg  = None
        if inrec[:8] != "KKoziaraEYWORDS":
            correct = False
            errmsg  = "KEYWORDS: first field should be KEYWORDS"
        else:
            keywords = [i.strip().lower() for i in inrec[8:].split(",")]
            for keyword in keywords:
                # stupid full stop
                if keyword[-1] == ".":
                    keyword = keyword[:-1]
                if keyword not in legal_keywords:
                    correct = False
                    errmsg  = "KEYWORDS: illegal keyword %s" % keyword
                    break
                else:
                    self.keywords.append(keyword)
        return correct, errmsg


    def _accuracy_check(self, inrec):
        '''
        Checking the accuracy field
        
        Input:
        inrec  :
        
        Output:
        [0]    : boolean   Correctness
        [1]    : String    Message
        '''
        
        correct = True
        errmsg  = None
        fields  = [i.strip() for i in inrec.split()]
        if len(fields) != 4:
            correct = False
            errmsg  = "ACCURACY: error in number of fields. Should be 4"
        elif fields[0] != "ACCURACY":
            correct = False
            errmsg  = "ACCURACY: first field should be 'ACCURACY'"
        elif not fields[1].isdigit() or len(fields[1]) != 1:
            correct = False
            errmsg  = "ACCURACY: second field should be a single digit"
        elif not pr_field.match(fields[2]):
            correct = False
            errmsg  = "ACCURACY: error in PR field"
        elif not rc_field.match(fields[3]):
            correct = False
            errmsg  = "ACCURACY: error in RC field"
        return correct, errmsg


    def _go_prediction_check(self, inrec):
        '''
        Checking the go_prediction field
        
        Input:
        inrec  :
        
        Output:
        [0]    : boolean   Correctness
        [1]    : String    Message
        '''
        
        correct = True
        errmsg  = None
        fields  = [i.strip() for i in inrec.split()]
        if len(fields) != 3:
            correct = False
            errmsg  = "GO prediction: wrong number of fields. Should be 3"
        elif not target_field.match(fields[0]):
            correct = False
            errmsg  = "GO prediction: error in first (Target ID) field"
        elif not go_field.match(fields[1]):
            correct = False
            errmsg  = "GO prediction: error in second (GO ID) field"
        elif not confidence_field.match(fields[2]):
            correct = False
            errmsg  = "GO prediction: error in third (confidence) field"
        elif float(fields[2]) > 1.0:
            correct = False
            errmsg  = "GO prediction: error in third (confidence) field. Cannot be > 1.0"
        else:
            self.data[fields[0]].append(
                {'term': fields[1], 'confidence': float(fields[2])})
        return correct, errmsg


    def _end_check(self, inrec):
        '''
        Checking the end field
        
        Input:
        inrec  :
        
        Output:
        [0]    : boolean   Correctness
        [1]    : String    Message
        '''
        
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 1:
            correct = False
            errmsg = "END: wrong number of fields. Should be 1"
        elif fields[0] != "END":
            correct = False
            errmsg = "END: record should include the word END only"
        return correct, errmsg


    def _handle_error(self, correct, errmsg, inrec):
        '''
        Error handler for predictions
        
        Input:
        correct : Boolean    Correctness
        errmsg  : String     Message
        inrec   : 
        '''
        
        if not correct:
            print(inrec)
            raise ValueError(errmsg)
            

    def read(self, pred_path):
        '''
        Read in the prediction
        
        Input:
        pred_path: String    Location of the prediction
        '''
        
        # Initilize Variables
        visited_states   = []
        n_accuracy       = 0
        first_prediction = True
        first_accuracy   = True
        first_keywords   = True
        n_models         = 0
        filename         = pred_path.name.split('/')[-1]
        filenamefields   = filename.split('.')[0].split('_')
        self.taxon       = filenamefields[2]
        #inline = open(pred_path)
        for inline in pred_path:
            # gzipped files are in bytes. Need to convert to utf-8
            if type(inline) is bytes:
                inline = inline.decode("utf-8")
            inrec = [i.strip() for i in inline.split()]
            field1 = inrec[0]
            # Check which field type (state) we are in
            if field1 == "AUTHOR":
                state = "author"
            elif field1 == "MODEL":
                state = "model"
            elif field1 == "KEYWORDS":
                state = "keywords"
            elif field1 == "ACCURACY":
                state = "accuracy"
            elif field1 == "END":
                state = "end"
            else:  # default to prediction state
                state = "go_prediction"
            # Check for errors according to state
            if state == "author":
                correct, errmsg = self._author_check(inline)
                self._handle_error(correct, errmsg, inline)
                visited_states.append(state)

            elif state == "model":
                n_models += 1
                n_accuracy = 0
                if n_models > 3:
                    raise ValueError("Too many models. Only up to 3 allowed")
                correct, errmsg = self._model_check(inline)
                if correct and self.model != int(filenamefields[1]):
                    correct = False
                    errmsg = 'MODEL: model number in file different from filename: %' % pred_path.name
                self._handle_error(correct, errmsg, inline)
                if n_models == 1:
                    visited_states.append(state)
            elif state == "keywords":
                if first_keywords:
                    visited_states.append(state)
                    first_keywords = False
                correct, errmsg = self._keywords_check(inline)
                self._handle_error(correct, errmsg, inline)
            elif state == "accuracy":
                if first_accuracy:
                    visited_states.append(state)
                    first_accuracy = False
                n_accuracy += 1
                if n_accuracy > 3:
                    self._handle_error(False, "ACCURACY: too many ACCURACY records: %" % pred_path.name)
                else:
                    correct, errmsg = self._accuracy_check(inline)
                    # Don't handle this error if its raised #################################################
            elif state == "go_prediction":
                correct, errmsg = self._go_prediction_check(inline)
                self._handle_error(correct, errmsg, inline)
                if first_prediction:
                    visited_states.append(state)
                    first_prediction = False
            elif state == "end":
                correct, errmsg = self._end_check(inline)
                self._handle_error(correct, errmsg, inline)
                visited_states.append(state)
            # End file forloop
        if (visited_states != legal_states1 and
            visited_states != legal_states2 and
            visited_states != legal_states3 and
            visited_states != legal_states4):
            print(visited_states)
            print("file not formatted according to CAFA specs")
            print("Check whether all these record types are in your file")
            print(
                "Check whether all these record types are in your file in the correct order")
            print("AUTHOR, MODEL, KEYWORDS, ACCURACY (optional), predictions, END")
            raise ValueError
            

    def read_and_split_and_write(self, obo_path, prediction_path):
        '''
        Description
        
        Input:
        obo_path        : file location of OBO file     
        prediction_path : should be a handle
        '''
        # This function has self.read() included
        # pred_path 
        # split both writes to the predictions to three separate files
        self.read(prediction_path)
        go_graph = OboIO.OboReader(open(obo_path)).read()
        mfo_out = open("%s_MFO.txt" % os.path.splitext(prediction_path.name)[0], "w")
        bpo_out = open("%s_BPO.txt" % os.path.splitext(prediction_path.name)[0], "w")
        cco_out = open("%s_CCO.txt" % os.path.splitext(prediction_path.name)[0], "w")
        for protein in self.data.items():
            for u in protein[1]:
                try:
                    namespace = go_graph.get_namespace(u['term'])
                except KeyError:
                    #print("Term %s not included in obo file.\n" % str(e))
                    namespace = None
                    # This excludes a term that's not found in the obo file
                    # to be entered into any of the ontology-specific files
                if namespace == 'molecular_function':
                    mfo_out.write("%s\t%s\t%.2f\n" %
                                  (protein[0], u['term'], u['confidence']))
                elif namespace == 'biological_process':
                    bpo_out.write("%s\t%s\t%.2f\n" %
                                  (protein[0], u['term'], u['confidence']))
                elif namespace == 'cellular_component':
                    cco_out.write("%s\t%s\t%.2f\n" %
                                  (protein[0], u['term'], u['confidence']))
                elif namespace == None:
                    #print("Some terms are not included in obo file\n")
                    continue
                else:
                    raise ValueError(
                        "Term %s not found in any ontology" % u['term'])
        mfo_out.close()
        bpo_out.close()
        cco_out.close()
