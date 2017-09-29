#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 13:04:14 2017

@author: pianarol
"""

import ldf
import numpy as np
import math

t = ldf.t
s = ldf.s

######################### For the time t

#make sure the array contains only float
t = np.array(ldf.t)
t = np.asfarray(t, float)


#take only the positive values
inds_t = [i for (i, val) in enumerate(t) if val > 0]
t = [val for (i, val) in enumerate(t) if i in inds_t]



######################### For the time t

#make sure the array contains only float
s = np.array(ldf.s)
s = np.asfarray(s, float)


#take only the positive values
inds_t = [i for (i, val) in enumerate(s) if val > 0]
s = [val for (i, val) in enumerate(s) if i in inds_t]

#create the table t1,t2,t3 and s1,s2,s3

endt = len(t)
ends = len(s)

t1 = t[0:endt-2]
t2 = t[1:endt-1]
t3 = t[2:endt]

endt1 = len(t1)

s1 = s[0:ends-2]
s2 = s[1:ends-1]
s3 = s[2:ends]


######from the hytool of matlab, to know what is what ##################
#d1 = (log(t2./t1).*s3)./       (log(t3./t2).*log(t3./t1));
#d2 = (log(t3.*t1./t2.^2).*s2)./(log(t3./t2).*log(t2./t1));
#d3 = (log(t3./t2).*s1)./ (log(t2./t1).*log(t3./t1));
#     

#### d1 ####

#log(t2/t1) 
logt2t1 = []
for i in range(0, endt1):
    logt2t1.append(math.log(t2[i]/t1[i]))

#log(t2/t1)*s3 : 
D1_part1 = np.array(logt2t1) * np.array(s3)
 
#log(t3/t2)
logt3t2 = []
for i in range(0, endt1):
    logt3t2.append(math.log(t3[i]/t2[i]))

#log(t3/t1)
logt3t1 = []
for i in range(0, endt1):
    logt3t1.append(math.log(t3[i]/t2[i]))

#log(t3/t2)*log(t3/t1)
D1_part2 = np.array(logt3t2)*np.array(logt3t1)

d1 = D1_part1/D1_part2

#### d2 ####

#log(t3*t1/t2²)
logt3t1t2 = []
for i in range(0, endt1):
    logt3t1t2.append(math.log(t3[i]*t1[i]/t2[i]**2))

#logt3t1t2 * s2
D2_part1 = np.array(logt3t1t2) * np.array(s2)

#log(t3/t2)*log(t2/t1)
D2_part2 = np.array(logt3t2)*np.array(logt2t1)

d2 = D2_part1 / D2_part2

#### d3 ####

#logt3/t2 * s1
D3_part1 = np.array(logt3t2) * np.array(s1)

#log(t2/t1)*log(t3/t2)
D3_part2 = np.array(logt2t1)*np.array(logt3t2)

d3 = D3_part1 / D3_part2

#xd and yd

xd = t2
yd = d1+d2-d3
