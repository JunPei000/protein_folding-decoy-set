#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 14:36:48 2017

@author: peijun
"""
from __future__ import division
from function import name
import json
from os import listdir
from math import exp,log, sqrt, pi
def tor_calc_para(dicT):
    R = 1.987*10**(-3)
    T = 298
    fpt = open('path of torion data set','r')
    Rpt = open('path of torsion range data set','r')
    T_data = json.load(fpt)
    Range_data = {}
    Range_data = json.load(Rpt)
    dicT_energy = {}
    for key in dicT:
        dicT_energy[key]=[]
        A = []
        A = name(key)
        if 'HIE' in A[0]:
            A[0] = A[0].replace('HIE','HID')
        if 'HIE' in A[1]:
            A[1] = A[1].replace('HIE','HID')
        if 'HIS' in A[0]:
            A[0] = A[0].replace('HIS','HID')
        if 'HIS' in A[1]:
            A[1] = A[1].replace('HIS','HID')
        if 'HIP' in A[0]:
            A[0] = A[0].replace('HIP','HID')
        if 'HIP' in A[1]:
            A[1] = A[1].replace('HIP','HID')
        if 'OXT' in A[0]:
            A[0] = A[0].replace('OXT','O')
        if 'OXT' in A[1]:
            A[1] = A[1].replace('OXT','O')
        ind1 = A[0].index('-')
        ind2 = A[1].index('-')
        if len(A[0][ind1+1:]) ==4 and (A[0][ind1+1]=='C' or A[0][ind1+1]=='N'):
            A[0] = A[0][:ind1+1]+A[0][ind1+2:]
        if len(A[1][ind2+1:]) ==4 and (A[1][ind2+1]=='C' or A[1][ind2+1]=='N'):
            A[1] = A[1][:ind2+1]+A[1][ind2+2:]   
        key1 = A[1]+'__'+A[0]
        key2 = A[0] +'__'+A[1]
        b1 = 0; b2 = 0; b3 = 0; b4 = 0; b5 = 0; b6 = 0; b7 = 0; b8 = 0; b9 = 0; b10 = 0
        b11 = 0; b12 = 0; blen = 0
        if T_data.has_key(key1):
            atom_pair_name = key1
        elif T_data.has_key(key2): 
            atom_pair_name = key2
        elif not T_data.has_key(key1) and not T_data.has_key(key2):
            print'error! missing torsion!', key1, key2
	    continue
        if T_data[atom_pair_name].has_key('b10'):
            b1 = T_data[atom_pair_name]['b1']
            b2 = T_data[atom_pair_name]['b2']
            b3 = T_data[atom_pair_name]['b3']
            b4 = T_data[atom_pair_name]['b4']
            b5 = T_data[atom_pair_name]['b5']
            b6 = T_data[atom_pair_name]['b6']
            b7 = T_data[atom_pair_name]['b7']
            b8 = T_data[atom_pair_name]['b8']
            b9 = T_data[atom_pair_name]['b9']
            b10 = T_data[atom_pair_name]['b10']
            b11 = T_data[atom_pair_name]['b11']
            b12 = T_data[atom_pair_name]['b12']
            blen = 12
        elif not T_data[atom_pair_name].has_key('b10') and T_data[atom_pair_name].has_key('b7'):
            b1 = T_data[atom_pair_name]['b1']
            b2 = T_data[atom_pair_name]['b2']
            b3 = T_data[atom_pair_name]['b3']
            b4 = T_data[atom_pair_name]['b4']
            b5 = T_data[atom_pair_name]['b5']
            b6 = T_data[atom_pair_name]['b6']
            b7 = T_data[atom_pair_name]['b7']
            b8 = T_data[atom_pair_name]['b8']
            b9 = T_data[atom_pair_name]['b9']
            blen = 9
        elif not T_data[atom_pair_name].has_key('b7') and T_data[atom_pair_name].has_key('b4'):
            b1 = T_data[atom_pair_name]['b1']
            b2 = T_data[atom_pair_name]['b2']
            b3 = T_data[atom_pair_name]['b3']
            b4 = T_data[atom_pair_name]['b4']
            b5 = T_data[atom_pair_name]['b5']
            b6 = T_data[atom_pair_name]['b6']
            blen = 6
        elif not T_data[atom_pair_name].has_key('b4') and T_data[atom_pair_name].has_key('b1'):
            b1 = T_data[atom_pair_name]['b1']
            b2 = T_data[atom_pair_name]['b2']
            b3 = T_data[atom_pair_name]['b3']
            blen = 3
        s = Range_data[atom_pair_name]['start']
        e = Range_data[atom_pair_name]['end']
        for i in range(len(dicT[key])):
            dist = dicT[key][i]
            r = dist%0.005
            r = round(r,3)
            if r <=0.0025:
                dist = dist-r
            elif r > 0.0025:
                dist = dist+(0.005-r)
            dist = round(dist,3)
            max_t = 6
            min_t = 1
            dist_vec = []
            if (max_t - dist) >= 0.1 and (dist-min_t)>= 0.1:
                for j in range(-20,21):
                    dist_vec.append(float("{0:.3f}".format(0.005*j+dist)))
            elif (max_t - dist) < 0.1:
                for j in range(-40,1):
                    dist_vec.append(float("{0:.3f}".format(0.005*j+dist)))
            elif (dist - min_t) < 0.1:
                for j in range(0,41):
                    dist_vec.append(float("{0:.3f}".format(0.005*j+dist)))
                    
            dist_vec_energy = []
            if len(dist_vec) != 41:
                print 'error! wrong dist_vec!'
                break
            elif len(dist_vec) == 41:
                for j in range(len(dist_vec)):
                    if dist_vec[j] >= s and dist_vec[j] <= e:
                        tor_p = 0
                        if blen == 3:
                            tor_p = b3/(b2*sqrt(2*pi)) * exp(-(dist_vec[j]-b1)**2 /(2*(b2**2)))
                            dist_vec_energy.append(tor_p)
                        elif blen == 6:
                            tor_p = b3/(b2*sqrt(2*pi)) * exp(-(dist_vec[j]-b1)**2 /(2*(b2**2)))+b6/(b5*sqrt(2*pi)) * exp(-(dist_vec[j]-b4)**2 /(2*(b5**2)))
                            dist_vec_energy.append(tor_p)
                        elif blen == 9:
                            tor_p = b3/(b2*sqrt(2*pi)) * exp(-(dist_vec[j]-b1)**2 /(2*(b2**2)))+b6/(b5*sqrt(2*pi)) * exp(-(dist_vec[j]-b4)**2 /(2*(b5**2)))+b9/(b8*sqrt(2*pi)) * exp(-(dist_vec[j]-b7)**2 /(2*(b8**2)))
                            dist_vec_energy.append(tor_p)
                        elif blen == 12:
                            tor_p = b3/(b2*sqrt(2*pi)) * exp(-(dist_vec[j]-b1)**2 /(2*(b2**2)))+b6/(b5*sqrt(2*pi)) * exp(-(dist_vec[j]-b4)**2 /(2*(b5**2)))+b9/(b8*sqrt(2*pi)) * exp(-(dist_vec[j]-b7)**2 /(2*(b8**2)))+b12/(b11*sqrt(2*pi)) * exp(-(dist_vec[j]-b10)**2 /(2*(b11**2)))
                            dist_vec_energy.append(tor_p)
			if tor_p <= 10**(-10):
			    tor_p = 10**(-10)
                    elif dist_vec[j] < s or dist_vec[j] > e:
                        tor_p = 0
                        tor_p = 10**(-10)
                        dist_vec_energy.append(tor_p)

            sum_dist_energy = sum(dist_vec_energy)
            dicT_energy[key].append(sum_dist_energy)
    dest = {}
    for key in dicT_energy:
        sum1 = 0
        for i in range(len(dicT_energy[key])):
	    if dicT_energy[key][i] <= 10**(-30):
                dicT_energy[key][i] = 10**(-30)
            sum1 += log(dicT_energy[key][i])
        dest[key] = sum1 - len(dicT_energy[key])*log(41)
    return dest

