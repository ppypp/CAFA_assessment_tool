import math
import pickle as cp
import os
import argparse
import sys
import yaml
import ICHelper
import numpy as np
import scipy as sc
import time

def main():

    for ontology in ['BPO','CCO','MFO']:
        counter = 0
        data = dict()
        IC = []
        # Read in IC Values
        reader = open("/home/mcgerten/Downloads/IC-MATLAB/{}.txt".format(ontology),"r")
        for line in reader:
            IC = line.split(",")
        reader.close()
        reader = open("/home/mcgerten/Downloads/IC-MATLAB/{}Terms.txt".format(ontology),"r")
        # For every annotation, create an id and add to data
        for line in reader:
            #print(line)
            fields = line.split(",")
            GOTerm = fields[0]
            try:
                #print(IC[counter])
                data[GOTerm] = [0,IC[counter-1]]
            except:
                print("THERE WAS A PROBLEM")
                data[GOTerm] = [0,0]
            #print data[counter]
            if (counter % 1000 == 0):
                print ("I have done {} annotations".format(counter))
            counter += 1
     
    
        print ("I have done {} annotations".format(counter))
        cp.dump(data, open("ICdata/ia_{}-L.map".format(ontology),"wb"))
        convertToReadable(data, "{}".format(ontology))
        


        
def convertToReadable(IC, ontology):
    '''
    Convert Map to human readable Values to manually check accuracy
    '''
    data  = open("ICdata/IAmap_{}.txt".format(ontology), 'w')
    for term in IC:
        data.write('{}\t {}\t {}\n'.format(term, IC[term][0], IC[term][1]))
    data.close


if __name__ == '__main__':
   main()
   