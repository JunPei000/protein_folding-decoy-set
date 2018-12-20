#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 11:08:18 2017

@author: peijun
"""
from __future__ import division
from os import listdir
from nonb_calc import nonb
from torsion_calc_para import tor_calc_para
from protein_info_finder import protein_info
from function import delete_H
from function import name
import pandas as pd
from ref_des_build import ref_des_build
from change_key_to_standard import change_key_to_standard as ckts
direct1 = 'paht of the decoy set'
newfile1 = open(direct1+'descriptort.csv','w')
newfile2 = open(direct1+'descriptorn.csv','w')
lfiles = listdir(direct1)
[t,n,dictor, dicnon] = ref_des_build()
for pdbname in lfiles:
    if '.pdb' in pdbname:
        print(pdbname)
        lfile = direct1+pdbname
        p_dele_H = delete_H(lfile)
        lprotein = protein_info(p_dele_H)
        lprotein_tor = lprotein[0]
        lprotein_nonb = lprotein[1]
        lprotein_tor_energy = {}
        lprotein_nonb_energy = {}
        lprotein_tor_energy = tor_calc_para(lprotein_tor)
        lprotein_nonb_energy = nonb(lprotein_nonb)
        lprotein_tor_energy = ckts(dictor, lprotein_tor_energy)
        lprotein_nonb_energy = ckts(dicnon, lprotein_nonb_energy)
        t[pdbname] = t['Pair_name'].map(lprotein_tor_energy)
        n[pdbname] = n['Pair_name'].map(lprotein_nonb_energy)
t.to_csv(newfile1)
n.to_csv(newfile2)
newfile1.close()
newfile2.close()        
