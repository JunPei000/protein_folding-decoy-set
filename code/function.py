# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 14:36:59 2015

@author: peijun
"""
from __future__ import division
from math import pi,sqrt
def delete_H(pdbfile):
    if '.gz' in pdbfile:
        with gzip.open(pdbfile,'rb') as f:
            orig_file = f.readlines()
    else:
        orig_file = file.readlines(open(pdbfile,'r'))
    output_file = []
    for line in orig_file:
        newlist = []
        newlist = line.split()
        if len(newlist) != 0:
            if newlist[0] == 'ATOM' and newlist[2][0] != 'H':
                output_file.append(line)
    return output_file
def distribution(x):
    Bin = []; Num = []; M = [];
    i = 0.0
    while i <= 21.0:
        Bin.append(i)
        i += 0.05
    for i in range(len(Bin)):
        dind = str(Bin[i]).index('.')
        Bin[i] = float(str(Bin[i])[:dind+4])     
    if x == []:
        return [0.0], [0.0]
    elif x != []:
        minx = min(x)
        maxx = max(x)
        for i in range(len(Bin)):
            if minx >= Bin[i] and minx <= Bin[i+1]:
                a = i
            if maxx >= Bin[i] and maxx <= Bin[i+1]:
                b = (i+1)
        for i in range(a, b):
            Num.append([])
        for i in range(len(x)):
            for j in range(a, b):
                if x[i] >= Bin[j] and x[i] <= Bin[j+1]:
                    Num[j-a].append(x[i])
        M.append(0.0)
        for i in range(len(Num)):
            if Num[i] == []:
                M.append(0.0)
            if Num[i] != []:
                m = len(Num[i])
                M.append(m)
        return Bin[a:b+1], M               
        
        
def atmtyperatio(x,a):
    if x == []:
        return 0.0
    if x != []:
        Count = x.count(a)
        m = Count/(len(x))
        return m


def intensity(x):
    Bin = []; Num = []; M = [];
    i = 1.0
    while i <= 21.0:
        Bin.append(i)
        i += 0.05
    for i in range(len(Bin)):
        dind = str(Bin[i]).index('.')
        Bin[i] = float(str(Bin[i])[:dind+4])     
    if x == []:
        return [0.0], [0.0]
    elif x != []:
        minx = min(x)
        maxx = max(x)
        for i in range(len(Bin)):
            if minx >= Bin[i] and minx <= Bin[i+1]:
                a = i
            if maxx >= Bin[i] and maxx <= Bin[i+1]:
                b = (i+1)
        for i in range(a, b):
            Num.append([])
        for i in range(len(x)):
            for j in range(a, b):
                if x[i] >= Bin[j] and x[i] <= Bin[j+1]:
                    Num[j-a].append(x[i])
        M.append(0.0)
        for i in range(len(Num)):
            if Num[i] == []:
                M.append(0.0)
            if Num[i] != []:
                m = len(Num[i])/(4/3*pi*((Bin[i+1])**3-(Bin[i])**3))
                M.append(m)
        return Bin[a:b+1], M


      
        
def distance(x,y):
    dis = sqrt((x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2)
    dis = str(dis)
    m = dis.index('.')
    dis = dis[:m+4]
    dis = float(dis)
    return dis
    
def check(x):
    x = str(x)
    a = x.index('.')
    m = x[:a+4]
    n = x[a+4:]
    if m !='' and n != '':
        return [m, n]
    if m == '' and n != '':
        return [n]
    if m != '' and n == '':
        return [m]

def tordel(x):
    new = []
    new.append(x[0])
    new.append(x[1])
    for i in range(2,len(x)):
        if x[i] not in new:
            a = x.count(x[i])
            if a == 2:
                new.append(x[i])
            if a > 2:
                b = int(a/2)
                for j in range(0,b):
                    new.append(x)
    return new

def name(x):
    a = x.index('_')
    A = x[:a]
    B = x[a+2:]
    return [A, B]  

def normalize(x):
    sum = 0
    for i in range(len(x)):
        x[i] = float(x[i])
        sum += x[i]
    for i in range(len(x)):
        x[i] = x[i]/sum
    return x
        
def solequa(A,B,C):
    a = float(A)
    b = float(B)
    c = float(C)
    if a == 0:
        return 'error! a equals one!'
    elif a != 0:
        d = b**2-4*a*c
        if d < 0:
            return 'error! d <0!'
        elif d >= 0:
            e = 0
            e = sqrt(d)
            s1 = (-b+e)/(2*a)
            s2 = (-b-e)/(2*a)
            s = []
            if s1 >0:
                s.append(s1)
            if s2 >0:
                s.append(s2)
            return s

def R(x,y):
    if len(x) != len(y):
        print ('error! Should give equal length array!')
    elif len(x) == len(y):
        A = []
        B = []
        A2 = []
        B2 = []
        C = []
        meanx = sum(x)/len(x)
        meany = sum(y)/len(y)
        for i in range(len(x)):
            A.append(x[i]-meanx)
            B.append(y[i]-meany)
            C.append((x[i]-meanx)*(y[i]-meany))
            A2.append((x[i]-meanx)**2)
            B2.append((y[i]-meany)**2)
        r = (sum(C))**2/(sum(A2)*sum(B2))
        return r
        
        
        
        
        
        
        
        
    
    

    
