# -*- coding: utf-8 -*-
"""
A helper class for functions needed throughout the program 
"""

def init(v_choice):
    '''
    Setup the global variables needed

    Input:
    v_choice : String    whether to be verbose
    '''
    
    global v
    
    if(v_choice == "Y"):
        v = 1
    elif(v_choice == "N"):
        v = 0
    else:
        # Throw an error?
        pass


# Could make this a tiered system -> differing numbers represent different verbosity
def vprint(s):
    '''
    Print only if Verbose is chosen
    
    Input:
    s : The string to print
    '''
    if v == 1 :
        print(s)