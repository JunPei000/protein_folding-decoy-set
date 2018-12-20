#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 14:13:33 2018

@author: peijun
"""
import numpy as np
import pandas as pd
import json
def ref_des_build():
    torsion = {}
    t = {}
    fpt = open('path of torsion data base','r')
    torsion = json.load(fpt)
    for key in torsion:
        if not t.has_key(key):
            t[key] = 1.0
    descriptort = pd.DataFrame(t.items(),columns=['Pair_name', 'ref'])
    
    nonb_file = open('path of nonbond data base','r')
    N_data = {}
    N_data = json.load(nonb_file)
    n = {}
    for key in N_data:
        if not n.has_key(key):
            n[key] = 1.0
    descriptorn = pd.DataFrame(n.items(), columns=['Pair_name', 'ref'])
    return descriptort, descriptorn,t,n



def ref_des_unif_100():
    torsion = {}
    t = {}
    fpt = open('/mnt/home/peijun/Documents/mlproject/protein_folding/KECSA2/unif_100/torsion_para_100.json','r')
    torsion = json.load(fpt)
    for key in torsion:
        if not t.has_key(key):
            t[key] = 1.0
    descriptort = pd.DataFrame(t.items(),columns=['Pair_name', 'ref'])
    
    
    nonb_file = open('/mnt/home/peijun/Documents/mlproject/protein_folding/KECSA2/unif_100/nonb_para_100.json','r')
    N_data = {}
    N_data = json.load(nonb_file)
    n = {}
    for key in N_data:
        if not n.has_key(key):
            n[key] = 1.0
    descriptorn = pd.DataFrame(n.items(), columns=['Pair_name', 'ref'])
    return descriptort, descriptorn,t,n    


def ref_des_unif_500():
    torsion = {}
    t = {}
    fpt = open('/mnt/home/peijun/Documents/mlproject/protein_folding/KECSA2/unif_500/torsion_para_500.json','r')
    torsion = json.load(fpt)
    for key in torsion:
        if not t.has_key(key):
            t[key] = 1.0
    descriptort = pd.DataFrame(t.items(),columns=['Pair_name', 'ref'])


    nonb_file = open('/mnt/home/peijun/Documents/mlproject/protein_folding/KECSA2/unif_500/nonb_para_500.json','r')
    N_data = {}
    N_data = json.load(nonb_file)
    n = {}
    for key in N_data:
        if not n.has_key(key):
            n[key] = 1.0
    descriptorn = pd.DataFrame(n.items(), columns=['Pair_name', 'ref'])
    return descriptort, descriptorn,t,n

def ref_des_unif_all():
    torsion = {}
    t = {}
    fpt = open('/mnt/home/peijun/Documents/mlproject/protein_folding/KECSA2/unif_all/torsion_para_all.json','r')
    torsion = json.load(fpt)
    for key in torsion:
        if key not in t:
            t[key] = 1.0
    descriptort = pd.DataFrame(t.items(),columns=['Pair_name', 'ref'])


    nonb_file = open('/mnt/home/peijun/Documents/mlproject/protein_folding/KECSA2/unif_all/nonb_para_all.json','r')
    N_data = {}
    N_data = json.load(nonb_file)
    n = {}
    for key in N_data:
        if key not in n:
            n[key] = 1.0
    descriptorn = pd.DataFrame(n.items(), columns=['Pair_name', 'ref'])
    return descriptort, descriptorn,t,n

