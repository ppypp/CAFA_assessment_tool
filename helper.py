# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 15:50:49 2018

@author: mcgerten
"""

def vprint(s):
    global v
    if v == 1 :
        # No newline character
        print(s)