#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 13:23:59 2018

@author: peijun
"""
####################################################################
#  This code is mainly about how to generate the final             #
#  descriptor for protein_folding machine learning project.        # 
#  The final descriptor should contain 2001+14028                  #
#  elements, for 2001 part, it contains all of the difference      #
#  between log(native torsion atompair A_B) and                    #
#  log(native torsion atompair A_B), and for the 14028 part, it    #
#  contains the nonbonding interaction, the last part contains     #
#  solvation information of native and decoy proteins              #
####################################################################


import os
import pandas as pd
    
torsionpart = pd.read_csv('path of the file contains descriptort')
nonbondpart = pd.read_csv('path of the file contains descriptorn')
nativename = 'name of the native pdb file'
newdescriptort = pd.DataFrame()
newdescriptorn = pd.DataFrame()
torsionpart = torsionpart.fillna(0)
nonbondpart = nonbondpart.fillna(0)
for column in torsionpart.columns:
    if 'Pair_name' in column:
        newdescriptort[column] = torsionpart[column].map(lambda x: x+'t')
        for row in nonbondpart.index:
            if nonbondpart.loc[row, column] == 0:
                nonbondpart.loc[row, column] = 'sasa'
        newdescriptorn[column] = nonbondpart[column].map(lambda x: x+'n')
    elif 'ref' in column:
        continue
    else:
        if column != nativename and 'Unnamed' not in column:
            zero_y = column.replace('.pdb', '000')
            one_y = column.replace('.pdb','111')
            newdescriptort[zero_y] = torsionpart[column]- torsionpart[nativename]
            newdescriptort[one_y] = torsionpart[nativename] - torsionpart[column]
            newdescriptorn[zero_y] = nonbondpart[column]- nonbondpart[nativename]
            newdescriptorn[one_y] = nonbondpart[nativename] - nonbondpart[column]
for column in torsionpart.columns:
    if column not in nonbondpart.columns:
        print (column)
for column in nonbondpart.columns:
    if column not in torsionpart.columns:
        print (column)
newdescriptort = newdescriptort.fillna(0)
newdescriptorn = newdescriptorn.fillna(0)
fp1 = open('path of final file', 'w')
m  = pd.DataFrame()
m = newdescriptort.append(newdescriptorn, ignore_index=True)
lr = len(m.index)
for column in m.columns:
    if column[-3:] == '111':
        m.loc[lr, column] = 1
    if column[-3:] == '000':
        m.loc[lr, column] = 0
    if 'Pair_name' in column:
        m.loc[lr,column] = 'class'
m = pd.DataFrame.transpose(m)
print (len(m.columns))
m.to_csv(fp1)
fp1.close()

        
        
        
        
        
        
      
        
        
        
        
        
        
        
