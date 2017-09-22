#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 15:23:37 2017

@author: pianarol
"""

import ths_dmo








filename = ths_dmo.filename 
global t
global s
  




      
with open(filename) as f:
    
    lines = f.readlines()
    t = [line.split()[0] for line in lines]
    s = [line.split()[1] for line in lines]

   


