# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 14:22:50 2017

@author: peijun
"""
from __future__ import division
from function import name
import json
from math import exp,log
def nonb(Interf):
    R = 1.987*10**(-3)
    T = 298
    nonb_file = open('path of the nonbond data base','r')
    N_data = {}
    N_data = json.load(nonb_file)
    dicN_energy = {}
    for key in Interf:
        B = []
        B = name(key)
        if 'HIE' in B[0]:
            B[0] = B[0].replace('HIE','HID')
        if 'HIE' in B[1]:
            B[1] = B[1].replace('HIE','HID')
        if 'HIS' in B[0]:
            B[0] = B[0].replace('HIS','HID')
        if 'HIS' in B[1]:
            B[1] = B[1].replace('HIS','HID')
        if 'HIP' in B[0]:
            B[0] = B[0].replace('HIP','HID')
        if 'HIP' in B[1]:
            B[1] = B[1].replace('HIP','HID')
        if 'OXT' in B[0]:
            B[0] = B[0].replace('OXT','O')
        if 'OXT' in B[1]:
            B[1] = B[1].replace('OXT','O')
        ind1 = B[0].index('-')
        ind2 = B[1].index('-')
        if len(B[0][ind1+1:]) ==4 and (B[0][ind1+1]=='C' or B[0][ind1+1]=='N'):
            B[0] = B[0][:ind1+1]+B[0][ind1+2:]
        if len(B[1][ind2+1:]) ==4 and (B[1][ind2+1]=='C' or B[1][ind2+1]=='N'):
            B[1] = B[1][:ind2+1]+B[1][ind2+2:]
        key3 = B[1]+'__'+B[0]
        key4 = B[0]+'__'+B[1]
        dicN_energy[key] = []
        if N_data.has_key(key4):
            E1 = float(N_data[key4]['E1'])
            E2 = float(N_data[key4]['E2'])
            a = float(N_data[key4]['a'])
            b = float(N_data[key4]['b'])
            sigma = float(N_data[key4]['sigma'])
            a1 = float(N_data[key4]['repa'][0])
            b1 = float(N_data[key4]['repa'][1])
            c1 = float(N_data[key4]['repa'][2])
            r_dist = float(N_data[key4]['x1'])
        elif N_data.has_key(key3):
            E1 = float(N_data[key3]['E1'])
            E2 = float(N_data[key3]['E2'])
            a = float(N_data[key3]['a'])
            b = float(N_data[key3]['b'])
            sigma = float(N_data[key3]['sigma'])
            a1 = float(N_data[key3]['repa'][0])
            b1 = float(N_data[key3]['repa'][1])
            c1 = float(N_data[key3]['repa'][2])
            r_dist = float(N_data[key3]['x1'])
        elif not N_data.has_key(key4) and not N_data.has_key(key3):
            print ('error! Cannot find nonb atom pair in database!', key4,key3, B[0], B[1], B[1][ind2+1])
        if N_data.has_key(key3) or N_data.has_key(key4):
            for i in range(len(Interf[key])):
                r = 0
                dist_vec_n = []
                r = float("{0:.3f}".format(Interf[key][i]))
                if (20 - r) >= 0.5 and (r - 2) >= 0.5:
                    for j in range(-100,101):
                        dist_vec_n.append(float("{0:.3f}".format(0.005*j+r)))
                elif (20 - r) < 0.5:
                    for j in range(-200,1):
                        dist_vec_n.append(float("{0:.3f}".format(0.005*j+r)))
                elif (r - 2) < 0.5:
                    for j in range(0,201):
                        dist_vec_n.append(float("{0:.3f}".format(0.005*j+r)))
                dist_vec_n_energy = []
                if len(dist_vec_n) != 201:
                    print ('error! Wrong dist_vec_n!')
                    break
                elif len(dist_vec_n) == 201:
                    for j in range(len(dist_vec_n)):
                        t = (E1*((sigma/dist_vec_n[j])**a)-E2*((sigma/dist_vec_n[j])**b))/(-R*T)
                        prob = 0
                        if dist_vec_n[j] >=2 and dist_vec_n[j]<=r_dist:
                            prob = a1*(dist_vec_n[j]**2)+b1*dist_vec_n[j]+c1
                        elif dist_vec_n[j] > r_dist:
                            if t <=500:
                                prob = exp((E1*((sigma/dist_vec_n[j])**a)-E2*((sigma/dist_vec_n[j])**b))/(-R*T))
                            elif t > 500:
                                prob = exp(500)
                        elif dist_vec_n[j] < 2:
                            prob = 10**(-30)
                        if prob <= 10**(-30):
                            prob = 10**(-30)
                        dist_vec_n_energy.append(prob)
                if dist_vec_n_energy != []:                
                    sum_dist_n_energy = 0
                    sum_dist_n_energy = sum(dist_vec_n_energy)
                    dicN_energy[key].append(sum_dist_n_energy)
                else:
                    dicN_energy[key].append(0)
    desn = {}
    for key in dicN_energy:
        sum2 = 0
        for i in range(len(dicN_energy[key])):
            sum2 += log(dicN_energy[key][i])
        desn[key] = sum2 - len(dicN_energy[key])*log(201)
    return desn
