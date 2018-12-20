# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 13:01:20 2017

@author: peijun
"""
from __future__ import division
from function import distance
from function import check
def protein_info(pdb):
    direct2 = 'path of parm94-2.mol2'
    metal = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra', 'Sc', 'Y', 'Ti', 
             'Zr', 'Hf', 'V', 'Nb', 'Ta', 'Cr', 'Mo', 'W', 'Mn', 'Tc', 'Re', 'Fe', 'Ru', 'Os', 'Co', 
             'Rh', 'Ir', 'Ni', 'Pd', 'Pt', 'Cu', 'Ag', 'Au', 'Zn', 'Cd', 'Hg', 'Al', 'Ga', 'In', 'Tl', 
             'Ge', 'Sn', 'Pb', 'Sb', 'Bi', 'Po']
    refer = file.readlines(open(direct2,'r'))
    for i in range(len(refer)):
        if '@<TRIPOS>ATOM' in refer[i]:
            a = i
        if '@<TRIPOS>BOND' in refer[i]:
            b = i
        if '@<TRIPOS>SUBSTRUCTURE' in refer[i]:
            c = i
    atom = []; bond= []
    for i in range(a+1,b):
        newlist = refer[i].split()
        atom.append([newlist[0], newlist[1], newlist[7]])
    for i in range(len(atom)): 
        if atom[i][1] == 'HB3' and atom[i-2][1] != 'HB1' and atom[i-1][1] == 'HB2':
            atom[i][1] = 'HB1'
        if atom[i][1] == 'HG3' and atom[i-2][1] != 'HG1' and atom[i-1][1] == 'HG2':
            atom[i][1] = 'HG1'
        if atom[i][1] == 'HA3' and atom[i-2][1] != 'HA1' and atom[i-1][1] == 'HA2':
            atom[i][1] = 'HA1'
        if atom[i][1] == 'HD3' and atom[i-2][1] != 'HD1' and atom[i-1][1] == 'HD2':
            atom[i][1] = 'HD1'
        if atom[i][1] == 'HE3' and atom[i-2][1] != 'HE1' and atom[i-1][1] == 'HE2':
            atom[i][1] = 'HE1'
        if atom[i][1] == 'HG13' and atom[i-2][1] != 'HG11' and atom[i-1][1] == 'HG12':
            atom[i][1] = 'HG11'
    
    for i in range(b+1,c):
        newlist = refer[i].split()
        bond.append([newlist[0], newlist[1],newlist[2]])
        
    #myfile = file.readlines(open(pdb, 'r'))
    myfile = pdb
    atomcoor = []; metalcoor = []
    for i in range(len(myfile)):
        newlist = []
        newlist = myfile[i].split()
        if 'ATOM' in newlist[0]:
            atomnum = newlist[1]
            atomtype = newlist[2]
            residuetype = newlist[3]
            if len(newlist) == 12:
                x = float(newlist[6])
                y = float(newlist[7])
                z = float(newlist[8])
                chainnum = newlist[4]
                num = newlist[5]
            elif len(newlist) != 12:
                if len(newlist[4]) == 1:
                    chainnum = newlist[4]
                    num = newlist[5]
                    a = check(newlist[6])
                    b = check(newlist[7])
                    if len(a) == 2 and len(b) == 1:
                        x = float(a[0])
                        c = check(a[1])
                        if len(c) == 1:
                            y = float(a[1])
                            z = float(newlist[7])
                        if len(c) == 2:
                            y = float(c[0])
                            z = float(c[1])
                    if len(a) == 1 and len(b) == 2:
                        x = float(newlist[6])
                        y = float(b[0])
                        z = float(b[1])
                    if len(a) == 1 and len(b) == 1:
                        x = float(newlist[6])
                        y = float(newlist[7])
                        z = float(newlist[8])
                if len(newlist[4]) != 1:
                    chainnum = newlist[4][0]
                    num = newlist[4][1:]
                    a = check(newlist[5])
                    b = check(newlist[6])
                    if len(a) == 2 and len(b) == 1:
                        x = float(a[0])
                        y = float(a[1])
                        z = float(newlist[7])
                    if len(a) == 1 and len(b) == 2:
                        x = float(a[0])
                        y = float(b[0])
                        z = float(b[1])
                    if len(a) == 1 and len(b) == 1:
                        x = float(newlist[5])
                        y = float(newlist[6])
                        z = float(newlist[7])
                    
                    
                        
            newatmtype = atomtype+'-'+residuetype
            atomcoor.append([atomtype, residuetype, num, newatmtype, x, y, z,atomnum, chainnum])
        if 'ATOM' in newlist[0] or 'HETATM' in newlist[0]:
            for j in range(len(metal)):
                if metal[j] not in newlist[-1]:
                    continue
                if metal[j] in newlist[-1]:
                    if ('HETATM' in newlist[0] and len(newlist) == 11) or ('ATOM' in newlist[0] and len(newlist) == 12):
                        matomtype = newlist[2]
                        mx = float(newlist[6])
                        my = float(newlist[7])
                        mz = float(newlist[8])
                    if len(newlist) != 11 and 'HETATM' in newlist[0]:
                        if newlist[0] != 'HETATM':
                            matomtype = newlist[1]
                            e = check(newlist[4])
                            f = check(newlist[5])
                            if len(e) ==2 and len(f) == 1:
                                mx = float(e[0])
                                g = check(e[1])
                                if len(g) ==1:
                                    my = float(e[1])
                                    mz = float(f[0])
                                if len(g) ==2:
                                    my = float(g[0])
                                    mz = float(g[1])
                            if len(e) ==1 and len(f) == 2:
                                mx = float(e[0])
                                my = float(f[0])
                                mz = float(f[1])
                            if len(e) == 1 and len(f) == 1:
                                mx = float(e[0])
                                my = float(f[0])
                                mz = float(newlist[6])
                        if newlist[0] == 'HETATM':
                            matomtype = newlist[2]
                            e = check(newlist[5])
                            f = check(newlist[6])
                            if len(e) == 2 and len(f) ==1:
                                mx = float(e[0])
                                g = check(e[1])
                                if len(g) == 2:
                                    my = float(g[0])
                                    mz = float(g[1])
                                if len(g) == 1:
                                    my = float(e[1])
                                    mz = float(f[0])
                            if len(e) == 1 and len(f) == 2:
                                mx = float(e[0])
                                my = float(f[0])
                                mz = float(f[1])
                    metalcoor.append([matomtype, mx, my, mz])
        """put atomcoor into different categories due to different residues"""
    residue = [];residuename=[];resindex = []
    for i in range(len(atomcoor)-1):
        if atomcoor[i][0] == 'N':
            m = atomcoor[i][1]
            residue.append([]) 
            residuename.append(m)
            resindex.append(i)
            m = 0
        else:
            continue
    for i in range(len(resindex)):
        if i != len(resindex)-1:            
            for j in range(resindex[i], resindex[i+1]):
                residue[i].append([atomcoor[j][3], atomcoor[j][4], atomcoor[j][5], atomcoor[j][6], atomcoor[j][0], atomcoor[j][1], atomcoor[j][7],atomcoor[j][2],atomcoor[j][8]])
        if i == len(resindex)-1:
            for j in range(resindex[i],len(atomcoor)):
                residue[i].append([atomcoor[j][3], atomcoor[j][4], atomcoor[j][5], atomcoor[j][6], atomcoor[j][0], atomcoor[j][1], atomcoor[j][7],atomcoor[j][2],atomcoor[j][8]])
        """delete residue which in 4 A around a metal ion"""
        """Rindex contains residue's index which in 4A around a metal ion, Rfinal contains non-repeated residue-index"""
    Rindex = []; Rfinal = []
    for i in range(len(metalcoor)):
        for j in range(len(residue)):
            for k in range(len(residue[j])):
                dis12 = 0
                atom1 = []; atom2 = []
                atom1.append(float(residue[j][k][1]))
                atom1.append(float(residue[j][k][2]))
                atom1.append(float(residue[j][k][3]))
                atom2.append(float(metalcoor[i][1]))
                atom2.append(float(metalcoor[i][2]))
                atom2.append(float(metalcoor[i][3]))
                dis12 = distance(atom1,atom2)
                if dis12 <= 4.0:
                    Rindex.append(j)
                    break
    for i in range(len(Rindex)):
        if Rindex[i] not in Rfinal:
            Rfinal.append(Rindex[i])
    Rfinal.sort()
    for i in range(1,len(Rfinal)+1):
        m = Rfinal[-i]
        residue.remove(residue[m])
        """find torsion and non-bonding information for each atom in every residue"""
        """assign number to each atom in different residues"""
    for i in range(len(residue)):
        new1 = []; new2= []; new3 = []; new4 = []; new5 = [];new6 =[]
        for j in range(len(residue[i])):
            if residue[i][j][4] == 'HN1':
                residue[i][j][4] = 'H1'
                residue[i][j][0] = residue[i][j][0].replace('HN1', 'H1')
            if residue[i][j][4] == 'HN2':
                residue[i][j][4] = 'H2'
                residue[i][j][0] = residue[i][j][0].replace('HN2', 'H2')
            if residue[i][j][4] == 'HN3':
                residue[i][j][4] = 'H3'
                residue[i][j][0] = residue[i][j][0].replace('HN3', 'H3')
        for j in range(len(residue[i])):
            new1.append(residue[i][j][4])
        if 'H1' in new1 and 'H2' in new1 and 'H3' in new1:
            for j in range(len(residue[i])):
                residue[i][j][5] = 'N'+residue[i][j][5]
        for j in range(len(residue[i])):
            new2.append(residue[i][j][4])
        if 'OXT' in new2:
            for j in range(len(residue[i])):
                residue[i][j][5] = 'C'+residue[i][j][5]
        for j in range(len(residue[i])):
            if residue[i][j][5] == 'HIS':
                new3.append(residue[i][j][4])
        if new3 != []:
            for j in range(len(residue[i])):
                if 'HE2' not in new3:
                    residue[i][j][5] = 'HID'
                if 'HD1' not in new3:
                    residue[i][j][5] = 'HIE'
                if 'HE2'  in new3 and 'HD1' in new3:
                    residue[i][j][5] = 'HIP'
        for j in range(len(residue[i])):
            if residue[i][j][5] == 'NHIS':
                new4.append(residue[i][j][4])
        if new4 != []:
            for j in range(len(residue[i])):
                if 'HE2' not in new4:
                    residue[i][j][5] = 'NHID'
                if 'HD1' not in new4:
                    residue[i][j][5] = 'NHIE'
                if 'HE2'  in new4 and 'HD1' in new4:
                    residue[i][j][5] = 'NHIP'
        for j in range(len(residue[i])):
            if residue[i][j][5] == 'CHIS':
                new5.append(residue[i][j][4])
        if new5 != []:
            for j in range(len(residue[i])):
                if 'HE2' not in new5:
                    residue[i][j][5] = 'CHID'
                if 'HD1' not in new5:
                    residue[i][j][5] = 'CHIE'
                if 'HE2'  in new5 and 'HD1' in new5:
                    residue[i][j][5] = 'CHIP'    
    for i in range(len(residue)):
        for j in range(len(residue[i])):
            residue[i][j][0] = residue[i][j][4]+'-'+residue[i][j][5]
    for i in range(len(residue)):
        for j in range(len(residue[i])):
            for k in range(len(atom)):
                if residue[i][j][4] == atom[k][1] and residue[i][j][5] == atom[k][2]:
                    residue[i][j].append(atom[k][0])
                    m = atom[k][0]+'_'+residue[i][j][7]+'_'+residue[i][j][8]
                    residue[i][j].append(m)
                    m = 0
    for i in range(len(residue)):
        j = 0
        while j <= (len(residue[i])-1):
            if len(residue[i][j]) != 11:
                print (residue[i][j])
                residue[i].remove(residue[i][j])
                j -= 1
            j += 1
    tr = 0
    rl = len(residue) -1
    while tr <= rl:
        if residue[tr] == []:
            residue.remove(residue[tr])
            tr -= 1
            rl -= 1
        tr += 1
        """find out torsion information for each atom in different residues"""
        """1. get the complete connectivity for each protein, find out N and C in every residue"""
    NCconn = []; ncconn = []
    for i in range(len(residue)):
        NCconn.append([])
        for j in range(len(residue[i])):
            if residue[i][j][4] == 'N':
                NCconn[i].append(residue[i][j][4])
                NCconn[i].append(residue[i][j][9])
            if residue[i][j][4] == 'C':
                NCconn[i].append(residue[i][j][4])
                NCconn[i].append(residue[i][j][9])
        """arrange the lists in NCconn, like ['N', n1, 'C', n2]"""
    for i in range(len(NCconn)):
        if len(NCconn[i]) == 4:
            if NCconn[i][0] == 'N' and NCconn[i][2] == 'C':
                ncconn.append(NCconn[i])
            if NCconn[i][0] == 'C' and NCconn[i][2] == 'N':
                ncconn.append([])
                m = []
                m.append(NCconn[i][2])
                m.append(NCconn[i][3])
                m.append(NCconn[i][0])
                m.append(NCconn[i][1])
                ncconn.append(m)
                
    allatomnum = [];Bond = []; Angle = []; Torsion = []; Nonbonding= []; Torsioncoor = []
    Bond1 = []; Angle1 = []; Torsion1 = []; fNonbonding = []; Bond2 = []; fNonbonding = []
    for i in range(len(residue)):
        for j in range(len(residue[i])):
            allatomnum.append(residue[i][j][9])
    newbond = []
    for i in range(len(bond)):
        if (bond[i][1] in allatomnum) and (bond[i][2] in allatomnum):
            newbond.append(bond[i])
        if (bond[i][1] not in allatomnum) or (bond[i][2] not in allatomnum):
            continue
    bond = newbond
    newbond = []
    for i in range(len(ncconn)-1):
        bondAD = []
        m = 1228+i
        bondAD = [m, ncconn[i][3],ncconn[i+1][1]]
        bond.append(bondAD)
    atomnumre = []
    for i in range(len(residue)):
        atomnumre.append([])
        atomnumre[i].append(residue[i][0][7])
        for j in range(len(residue[i])):
            atomnumre[i].append(residue[i][j][9])
    for i in range(len(residue)):
        Bond2.append([])
        Angle.append([])
        Bond1.append([])
        Angle1.append([])
        Torsion.append([])
        Torsion1.append([])
        Nonbonding.append([])
        fNonbonding.append([])
        for j in range(len(residue[i])):
            Bond2[i].append([])
            Angle[i].append([])
            Bond1[i].append([])
            Angle1[i].append([])
            Torsion[i].append([])
            Torsion1[i].append([])
            Nonbonding[i].append([])
            fNonbonding[i].append([])
            for k in range(len(bond)):
                if residue[i][j][9] == bond[k][1] and bond[k][2] in atomnumre[i][1:]:
                    for l in range(len(residue[i])):
                        if residue[i][l][9] == bond[k][2]:
                            a = residue[i][l][1]
                            b = residue[i][l][2]
                            c = residue[i][l][3] 
                            d = residue[i][l][0]
                            e = residue[i][l][4]
                    if [residue[i][j][-1], bond[k][2]+'_'+residue[i][0][7]+'_'+residue[i][0][8], residue[i][j][1], residue[i][j][2], residue[i][j][3],a,b,c,residue[i][j][0],d] not in Bond2[i][j]:                     
                        if not ((residue[i][j][4] == 'C' and e == 'N') or (residue[i][j][4] == 'N' and e == 'C')):                             
                            Bond1[i][j].append(bond[k][2]+'_'+residue[i][0][7]+'_'+residue[i][0][8])
                            Bond2[i][j].append([residue[i][j][-1], bond[k][2]+'_'+residue[i][0][7]+'_'+residue[i][0][8], residue[i][j][1], residue[i][j][2], residue[i][j][3],a,b,c,residue[i][j][0],d])
                if residue[i][j][9] == bond[k][2] and bond[k][1] in atomnumre[i][1:]:
                    for l in range(len(residue[i])):
                        if residue[i][l][9] == bond[k][1]:
                            a = residue[i][l][1]
                            b = residue[i][l][2]
                            c = residue[i][l][3]
                            d = residue[i][l][0]
                            e = residue[i][l][4]
                    if [residue[i][j][-1], bond[k][1]+'_'+residue[i][0][7]+'_'+residue[i][0][8], residue[i][j][1], residue[i][j][2], residue[i][j][3],a,b,c,residue[i][j][0],d] not in Bond2[i][j]:                     
                        if not ((residue[i][j][4] == 'C' and e == 'N') or (residue[i][j][4] == 'N' and e == 'C')):                             
                            Bond2[i][j].append([residue[i][j][-1], bond[k][1]+'_'+residue[i][0][7]+'_'+residue[i][0][8], residue[i][j][1], residue[i][j][2], residue[i][j][3],a,b,c,residue[i][j][0],d])
                            Bond1[i][j].append(bond[k][1]+'_'+residue[i][0][7]+'_'+residue[i][0][8])
                if (len(residue[i][j][5]) == 3 and i<(len(residue) -1)) or (len(residue[i][j][5]) == 4 and residue[i][j][5][0] == 'N' and i<(len(residue) -1)):
                    s = 0
                    s = abs(int(atomnumre[i+1][0]) - int(atomnumre[i][0]))
                    if s == 1:     
                        if residue[i][j][9] == bond[k][1] and bond[k][2] in atomnumre[i+1][1:]:
                            for l in range(len(residue[i+1])):
                                if residue[i+1][l][9] == bond[k][2]:
                                    a = residue[i+1][l][1]
                                    b = residue[i+1][l][2]
                                    c = residue[i+1][l][3]
                                    d = residue[i+1][l][0]
                                    e = residue[i+1][l][4]
                            if [residue[i][j][-1], bond[k][2]+'_'+residue[i+1][0][7]+'_'+residue[i+1][0][8], residue[i][j][1], residue[i][j][2], residue[i][j][3],a,b,c,residue[i][j][0],d] not in Bond2[i][j]:  
                                if (residue[i][j][4] == 'C' and e == 'N'):                             
                                    Bond1[i][j].append(bond[k][2]+'_'+residue[i+1][0][7]+'_'+residue[i+1][0][8])
                                    Bond2[i][j].append([residue[i][j][-1], bond[k][2]+'_'+residue[i+1][0][7]+'_'+residue[i+1][0][8], residue[i][j][1], residue[i][j][2], residue[i][j][3],a,b,c,residue[i][j][0],d])
                        if residue[i][j][9] == bond[k][2] and bond[k][1] in atomnumre[i+1][1:]:
                            for l in range(len(residue[i+1])):
                                if residue[i+1][l][9] == bond[k][1]:
                                    a = residue[i+1][l][1]
                                    b = residue[i+1][l][2]
                                    c = residue[i+1][l][3]
                                    d = residue[i+1][l][0]
                                    e = residue[i+1][l][4]
                            if [residue[i][j][-1], bond[k][1]+'_'+residue[i+1][0][7]+'_'+residue[i+1][0][8], residue[i][j][1], residue[i][j][2], residue[i][j][3],a,b,c,residue[i][j][0],d] not in Bond2[i][j]:  
                                if (residue[i][j][4] == 'C' and e == 'N'):                                                             
                                    Bond2[i][j].append([residue[i][j][-1], bond[k][1]+'_'+residue[i+1][0][7]+'_'+residue[i+1][0][8], residue[i][j][1], residue[i][j][2], residue[i][j][3],a,b,c,residue[i][j][0],d])
                                    Bond1[i][j].append(bond[k][1]+'_'+residue[i+1][0][7]+'_'+residue[i+1][0][8])
                if (len(residue[i][j][5]) == 3 and i> 0) or (len(residue[i][j][5]) == 4 and residue[i][j][5][0] == 'C' and i> 0):
                    s = 0
                    s = abs(int(atomnumre[i-1][0]) - int(atomnumre[i][0]))
                    if s == 1:
                        if residue[i][j][9] == bond[k][1] and bond[k][2] in atomnumre[i-1][1:]:
                            for l in range(len(residue[i-1])):
                                if residue[i-1][l][9] == bond[k][2]:
                                    a = residue[i-1][l][1]
                                    b = residue[i-1][l][2]
                                    c = residue[i-1][l][3]
                                    d = residue[i-1][l][0]
                                    e = residue[i-1][l][4]
                            if [residue[i][j][-1], bond[k][2]+'_'+residue[i-1][0][7]+'_'+residue[i-1][0][8], residue[i][j][1], residue[i][j][2], residue[i][j][3],a,b,c,residue[i][j][0],d] not in Bond2[i][j]:  
                                if (residue[i][j][4] == 'N' and e == 'C'):                                                             
                                    Bond1[i][j].append(bond[k][2]+'_'+residue[i-1][0][7]+'_'+residue[i-1][0][8])
                                    Bond2[i][j].append([residue[i][j][-1], bond[k][2]+'_'+residue[i-1][0][7]+'_'+residue[i-1][0][8], residue[i][j][1], residue[i][j][2], residue[i][j][3],a,b,c,residue[i][j][0],d])
                        if residue[i][j][9] == bond[k][2] and bond[k][1] in atomnumre[i-1][1:]:
                            for l in range(len(residue[i-1])):
                                if residue[i-1][l][9] == bond[k][1]:
                                    a = residue[i-1][l][1]
                                    b = residue[i-1][l][2]
                                    c = residue[i-1][l][3]
                                    d = residue[i-1][l][0]
                                    e = residue[i-1][l][4]
                            if [residue[i][j][-1], bond[k][1]+'_'+residue[i-1][0][7]+'_'+residue[i-1][0][8], residue[i][j][1], residue[i][j][2], residue[i][j][3],a,b,c,residue[i][j][0],d] not in Bond2[i][j]:  
                                if (residue[i][j][4] == 'N' and e == 'C'):                                                                                             
                                    Bond2[i][j].append([residue[i][j][-1], bond[k][1]+'_'+residue[i-1][0][7]+'_'+residue[i-1][0][8], residue[i][j][1], residue[i][j][2], residue[i][j][3],a,b,c,residue[i][j][0],d])
                                    Bond1[i][j].append(bond[k][1]+'_'+residue[i-1][0][7]+'_'+residue[i-1][0][8])
                
    for i in range(len(Bond2)):
        for j in range(len(Bond2[i])):
            for k in range(len(Bond2[i][j])):
                if len(Bond2[i][j][k]) == 10:
                    for l in range(len(Bond2[i])):
                        for m in range(len(Bond2[i][l])):
                            if Bond2[i][j][k][1] == Bond2[i][l][m][0] and Bond2[i][l][m][1] not in Bond2[i][j][k]:
                                if [Bond2[i][j][k][0], Bond2[i][j][k][1], Bond2[i][l][m][1], Bond2[i][j][k][2], Bond2[i][j][k][3], Bond2[i][j][k][4]] not in Angle[i][j]:
                                    Angle[i][j].append([Bond2[i][j][k][0], Bond2[i][j][k][1], Bond2[i][l][m][1], Bond2[i][j][k][2], Bond2[i][j][k][3], Bond2[i][j][k][4]])
                                    Angle1[i][j].append(Bond2[i][l][m ][1])
                            if Bond2[i][j][k][1] == Bond2[i][l][m][1] and Bond2[i][l][m][0] not in Bond2[i][j][k]:
                                if [Bond2[i][j][k][0], Bond2[i][j][k][1], Bond2[i][l][m][0], Bond2[i][j][k][2], Bond2[i][j][k][3], Bond2[i][j][k][4]] not in Angle[i][j]:                                    
                                    Angle[i][j].append([Bond2[i][j][k][0], Bond2[i][j][k][1], Bond2[i][l][m][0], Bond2[i][j][k][2], Bond2[i][j][k][3], Bond2[i][j][k][4]]) 
                                    Angle1[i][j].append(Bond2[i][l][m][0])
                if (len(Bond2[i][j][k]) == 10 and len(residue[i][j][5]) == 3 and i<(len(residue) -1)) or (len(Bond2[i][j][k]) == 10 and len(residue[i][j][5]) == 4 and residue[i][j][5][0] == 'N' and i<(len(residue) -1)):
                    s = 0
                    s = abs(int(atomnumre[i+1][0]) - int(atomnumre[i][0]))
                    if s == 1:  
                        for l in range(len(Bond2[i+1])):
                            for m in range(len(Bond2[i+1][l])):
                                if Bond2[i][j][k][1] == Bond2[i+1][l][m][0] and Bond2[i+1][l][m][1] not in Bond2[i][j][k]:
                                    if [Bond2[i][j][k][0], Bond2[i][j][k][1], Bond2[i+1][l][m][1], Bond2[i][j][k][2], Bond2[i][j][k][3], Bond2[i][j][k][4]] not in Angle[i][j]:
                                        Angle[i][j].append([Bond2[i][j][k][0], Bond2[i][j][k][1], Bond2[i+1][l][m][1], Bond2[i][j][k][2], Bond2[i][j][k][3], Bond2[i][j][k][4]])
                                        Angle1[i][j].append(Bond2[i+1][l][m ][1])
                                if Bond2[i][j][k][1] == Bond2[i+1][l][m][1] and Bond2[i+1][l][m][0] not in Bond2[i][j][k]:
                                    if [Bond2[i][j][k][0], Bond2[i][j][k][1], Bond2[i+1][l][m][0], Bond2[i][j][k][2], Bond2[i][j][k][3], Bond2[i][j][k][4]] not in Angle[i][j]:                                    
                                        Angle[i][j].append([Bond2[i][j][k][0], Bond2[i][j][k][1], Bond2[i+1][l][m][0], Bond2[i][j][k][2], Bond2[i][j][k][3], Bond2[i][j][k][4]]) 
                                        Angle1[i][j].append(Bond2[i+1][l][m][0])
                if (len(Bond2[i][j][k]) == 10 and len(residue[i][j][5]) == 3 and i> 0) or (len(Bond2[i][j][k]) == 10 and len(residue[i][j][5]) == 4 and residue[i][j][5][0] == 'C' and i> 0):
                    s = 0
                    s = abs(int(atomnumre[i-1][0]) - int(atomnumre[i][0]))
                    if s == 1:
                        for l in range(len(Bond2[i-1])):
                            for m in range(len(Bond2[i-1][l])):
                                if Bond2[i][j][k][1] == Bond2[i-1][l][m][0] and Bond2[i-1][l][m][1] not in Bond2[i][j][k]:
                                    if [Bond2[i][j][k][0], Bond2[i][j][k][1], Bond2[i-1][l][m][1], Bond2[i][j][k][2], Bond2[i][j][k][3], Bond2[i][j][k][4]] not in Angle[i][j]:
                                        Angle[i][j].append([Bond2[i][j][k][0], Bond2[i][j][k][1], Bond2[i-1][l][m][1], Bond2[i][j][k][2], Bond2[i][j][k][3], Bond2[i][j][k][4]])
                                        Angle1[i][j].append(Bond2[i-1][l][m ][1])
                                if Bond2[i][j][k][1] == Bond2[i-1][l][m][1] and Bond2[i-1][l][m][0] not in Bond2[i][j][k]:
                                    if [Bond2[i][j][k][0], Bond2[i][j][k][1], Bond2[i-1][l][m][0], Bond2[i][j][k][2], Bond2[i][j][k][3], Bond2[i][j][k][4]] not in Angle[i][j]:
                                        Angle[i][j].append([Bond2[i][j][k][0], Bond2[i][j][k][1], Bond2[i-1][l][m][0], Bond2[i][j][k][2], Bond2[i][j][k][3], Bond2[i][j][k][4]]) 
                                        Angle1[i][j].append(Bond2[i-1][l][m][0])
    for i in range(len(Angle)):
        for j in range(len(Angle[i])):
            for k in range(len(Angle[i][j])):
                if len(Angle[i][j][k]) == 6:
                    for l in range(len(Bond2[i])):
                        for m in range(len(Bond2[i][l])):
                            if Angle[i][j][k][2] == Bond2[i][l][m][0] and Bond2[i][l][m][1] not in Angle[i][j][k]:
                                if [Angle[i][j][k][0], Angle[i][j][k][1], Angle[i][j][k][2], Bond2[i][l][m][1], Angle[i][j][k][3], Angle[i][j][k][4], Angle[i][j][k][5], Bond2[i][l][m][5], Bond2[i][l][m][6], Bond2[i][l][m][7], Bond2[i][l][m][9]] not in Torsion[i][j]:
                                    Torsion[i][j].append([Angle[i][j][k][0], Angle[i][j][k][1], Angle[i][j][k][2], Bond2[i][l][m][1], Angle[i][j][k][3], Angle[i][j][k][4], Angle[i][j][k][5], Bond2[i][l][m][5], Bond2[i][l][m][6], Bond2[i][l][m][7], Bond2[i][l][m][9]])
                                    Torsion1[i][j].append(Bond2[i][l][m ][1])
                            if Angle[i][j][k][2] == Bond2[i][l][m][1] and Bond2[i][l][m][0] not in Angle[i][j][k]:
                                if [Angle[i][j][k][0], Angle[i][j][k][1], Angle[i][j][k][2], Bond2[i][l][m][0], Angle[i][j][k][3], Angle[i][j][k][4], Angle[i][j][k][5], Bond2[i][l][m][2], Bond2[i][l][m][3], Bond2[i][l][m][4], Bond2[i][l][m][8]] not in Torsion[i][j]:                                        
                                    Torsion[i][j].append([Angle[i][j][k][0], Angle[i][j][k][1], Angle[i][j][k][2], Bond2[i][l][m][0], Angle[i][j][k][3], Angle[i][j][k][4], Angle[i][j][k][5], Bond2[i][l][m][2], Bond2[i][l][m][3], Bond2[i][l][m][4], Bond2[i][l][m][8]]) 
                                    Torsion1[i][j].append(Bond2[i][l][m][0])
                if (len(Angle[i][j][k]) == 6 and len(residue[i][j][5]) == 3 and i<(len(residue) -1)) or (len(Angle[i][j][k]) == 6 and len(residue[i][j][5]) == 4 and residue[i][j][5][0] == 'N' and i<(len(residue) -1)):
                    s = 0
                    s = abs(int(atomnumre[i+1][0]) - int(atomnumre[i][0]))
                    if s == 1:  
                        for l in range(len(Bond2[i+1])):
                            for m in range(len(Bond2[i+1][l])):
                                if Angle[i][j][k][2] == Bond2[i+1][l][m][0] and Bond2[i+1][l][m][1] not in Angle[i][j][k]:
                                    if [Angle[i][j][k][0], Angle[i][j][k][1], Angle[i][j][k][2], Bond2[i+1][l][m][1], Angle[i][j][k][3], Angle[i][j][k][4], Angle[i][j][k][5], Bond2[i+1][l][m][5], Bond2[i+1][l][m][6], Bond2[i+1][l][m][7], Bond2[i+1][l][m][9]] not in Torsion[i][j]:
                                        Torsion[i][j].append([Angle[i][j][k][0], Angle[i][j][k][1], Angle[i][j][k][2], Bond2[i+1][l][m][1], Angle[i][j][k][3], Angle[i][j][k][4], Angle[i][j][k][5], Bond2[i+1][l][m][5], Bond2[i+1][l][m][6], Bond2[i+1][l][m][7], Bond2[i+1][l][m][9]])
                                        Torsion1[i][j].append(Bond2[i+1][l][m][1])
                                if Angle[i][j][k][2] == Bond2[i+1][l][m][1] and Bond2[i+1][l][m][0] not in Angle[i][j][k]:
                                    if [Angle[i][j][k][0], Angle[i][j][k][1], Angle[i][j][k][2], Bond2[i+1][l][m][0], Angle[i][j][k][3], Angle[i][j][k][4], Angle[i][j][k][5], Bond2[i+1][l][m][2], Bond2[i+1][l][m][3], Bond2[i+1][l][m][4], Bond2[i+1][l][m][8]] not in Torsion[i][j]:                                    
                                        Torsion[i][j].append([Angle[i][j][k][0], Angle[i][j][k][1], Angle[i][j][k][2], Bond2[i+1][l][m][0], Angle[i][j][k][3], Angle[i][j][k][4], Angle[i][j][k][5], Bond2[i+1][l][m][2], Bond2[i+1][l][m][3], Bond2[i+1][l][m][4], Bond2[i+1][l][m][8]]) 
                                        Torsion1[i][j].append(Bond2[i+1][l][m][0])
                if (len(Angle[i][j][k]) == 6 and len(residue[i][j][5]) == 3 and i> 0) or (len(Angle[i][j][k]) == 6 and len(residue[i][j][5]) == 4 and residue[i][j][5][0] == 'C' and i> 0):
                    s = 0
                    s = abs(int(atomnumre[i-1][0]) - int(atomnumre[i][0]))
                    if s == 1:
                        for l in range(len(Bond2[i-1])):
                            for m in range(len(Bond2[i-1][l])):
                                if Angle[i][j][k][2] == Bond2[i-1][l][m][0] and Bond2[i-1][l][m][1] not in Angle[i][j][k]:
                                    if [Angle[i][j][k][0], Angle[i][j][k][1], Angle[i][j][k][2], Bond2[i-1][l][m][1], Angle[i][j][k][3], Angle[i][j][k][4], Angle[i][j][k][5], Bond2[i-1][l][m][5], Bond2[i-1][l][m][6], Bond2[i-1][l][m][7], Bond2[i-1][l][m][9]] not in Torsion[i][j]:                                    
                                        Torsion[i][j].append([Angle[i][j][k][0], Angle[i][j][k][1], Angle[i][j][k][2], Bond2[i-1][l][m][1], Angle[i][j][k][3], Angle[i][j][k][4], Angle[i][j][k][5], Bond2[i-1][l][m][5], Bond2[i-1][l][m][6], Bond2[i-1][l][m][7], Bond2[i-1][l][m][9]])
                                        Torsion1[i][j].append(Bond2[i-1][l][m ][1])
                                if Angle[i][j][k][2] == Bond2[i-1][l][m][1] and Bond2[i-1][l][m][0] not in Angle[i][j][k]:
                                    if [Angle[i][j][k][0], Angle[i][j][k][1], Angle[i][j][k][2], Bond2[i-1][l][m][0], Angle[i][j][k][3], Angle[i][j][k][4], Angle[i][j][k][5], Bond2[i-1][l][m][2], Bond2[i-1][l][m][3], Bond2[i-1][l][m][4], Bond2[i-1][l][m][8]] not in Torsion[i][j]:                                    
                                        Torsion[i][j].append([Angle[i][j][k][0], Angle[i][j][k][1], Angle[i][j][k][2], Bond2[i-1][l][m][0], Angle[i][j][k][3], Angle[i][j][k][4], Angle[i][j][k][5], Bond2[i-1][l][m][2], Bond2[i-1][l][m][3], Bond2[i-1][l][m][4], Bond2[i-1][l][m][8]]) 
                                        Torsion1[i][j].append(Bond2[i-1][l][m][0]) 
                                        
    finalTatmtype = [];finalTatmnum = []; finalTlength = []
    for i in range(len(residue)):
        for j in range(len(residue[i])):
            if residue[i][j][0] not in finalTatmtype:
                finalTatmtype.append(residue[i][j][0])
    for i in range(len(finalTatmtype)):
        finalTatmnum.append([])
        finalTlength.append([])
        for j in range(len(residue)):
            for k in range(len(residue[j])):
                if residue[j][k][0] == finalTatmtype[i]:
                    m = residue[j][k][-1]
                    if m not in finalTatmnum[i]:
                        finalTatmnum[i].append(m)
                    m = 0
    for i in range(len(Torsion)):
        for j in range(len(Torsion[i])):
            for k in range(len(finalTatmnum)):
                if Torsion[i][j] != []:
                    if Torsion[i][j][0][0] in finalTatmnum[k]:
                        finalTlength[k].append(Torsion[i][j])
    for i in range(len(finalTlength)):
        for j in range(len(finalTlength[i])):
            for k in range(len(finalTlength[i][j])):
                X =[]; Y = []
                X.append(float(finalTlength[i][j][k][4]))
                X.append(float(finalTlength[i][j][k][5]))
                X.append(float(finalTlength[i][j][k][6]))
                Atype = finalTatmtype[i]
                Y.append(float(finalTlength[i][j][k][7]))
                Y.append(float(finalTlength[i][j][k][8]))
                Y.append(float(finalTlength[i][j][k][9]))
                Btype = finalTlength[i][j][k][10]
                dist = distance(X,Y)
                finalTlength[i][j][k].append(Atype)
                finalTlength[i][j][k].append(Btype)
                finalTlength[i][j][k].append(dist)
    finalT = []
    for i in range(len(finalTlength)):
        finalT.append([])
        new = []
        for j in range(len(finalTlength[i])):
            for k in range(len(finalTlength[i][j])):
                if finalTlength[i][j][k][-2] not in new:
                    new.append(finalTlength[i][j][k][-2])
        for k in range(len(new)):
            finalT[i].append([])
            finalT[i][k].append(finalTlength[i][j][0][-3])
            finalT[i][k].append(new[k])
        for j in range(len(finalTlength[i])):
            for k in range(len(finalTlength[i][j])):
                for l in range(len(finalT[i])):
                    if finalTlength[i][j][k][-3] == finalT[i][l][0] and finalTlength[i][j][k][-2] == finalT[i][l][1]:
                        finalT[i][l].append(finalTlength[i][j][k][-1])
    for i in range(len(residue)):
        for j in range(len(residue[i])):
            x1 = float(residue[i][j][1])
            y1 = float(residue[i][j][2])
            z1 = float(residue[i][j][3])
            atmn1 = residue[i][j][-1]
            atmt1 = residue[i][j][0]
            A = [x1,y1,z1]
            for k in range(len(residue)):
                for l in range(len(residue[k])):
                    if (k == i and l >= j) or k > i:
                        x2 = float(residue[k][l][1])
                        y2 = float(residue[k][l][2])
                        z2 = float(residue[k][l][3])
                        atmn2 = residue[k][l][-1]
                        atmt2 = residue[k][l][0]
                        B = [x2, y2, z2]                        
                        if abs(x2-x1) <= 20 and abs(y2-y1) <= 20 and abs(z2-z1) <= 20:
                            dis = distance(A,B)
                            if dis <= 30:
                                Nonbonding[i][j].append([atmt1,atmn1,atmt2,atmn2,dis])
    for i in range(len(Nonbonding)):
        for j in range(len(Nonbonding[i])):
            for k in range(len(Nonbonding[i][j])):
                if Nonbonding[i][j][k][-1] != 0 and Nonbonding[i][j][k][0][0] != 'H' and Nonbonding[i][j][k][2][0] != 'H':
                    if Nonbonding[i][j][k][3] not in Bond1[i][j] and Nonbonding[i][j][k][3] not in Angle1[i][j] and Nonbonding[i][j][k][3] not in Torsion1[i][j]:
                        fNonbonding[i][j].append(Nonbonding[i][j][k])
    dicT = {}
    dicN = {}
    dicN_info = {}
    for i in range(len(finalT)):
        for j in range(len(finalT[i])):
            if finalT[i][j][0][0] != 'H' and finalT[i][j][1][0] != 'H' and finalT[i][j][0] != finalT[i][j][1]:
                atomp_name1 = finalT[i][j][0]+'__'+finalT[i][j][1]
                atomp_name2 = finalT[i][j][1]+'__'+finalT[i][j][0]
                if dicT.has_key(atomp_name1):
                    for k in range(2,len(finalT[i][j])):
                        dicT[atomp_name1].append(finalT[i][j][k])
                elif not dicT.has_key(atomp_name1) and dicT.has_key(atomp_name2):
                    continue
                elif not dicT.has_key(atomp_name1) and not dicT.has_key(atomp_name2):
                    dicT[atomp_name1] = []
                    for k in range(2,len(finalT[i][j])):
                        dicT[atomp_name1].append(finalT[i][j][k])
            elif finalT[i][j][0][0] != 'H' and finalT[i][j][1][0] != 'H' and finalT[i][j][0] == finalT[i][j][1]:
                atomp_name = finalT[i][j][0]+'__'+finalT[i][j][1]
                if dicT.has_key(atomp_name):
                    print ('error!', atomp_name)
                elif not dicT.has_key(atomp_name):
                    dicT[atomp_name] = []
                    for k in range(2,len(finalT[i][j])):
                        if finalT[i][j][k] not in dicT[atomp_name]:
                            dicT[atomp_name].append(finalT[i][j][k])
                        else:
                            continue
                    if len(dicT[atomp_name])!=((len(finalT[i][j])-2)/2):
                        #print 'error!, did not paste enough data!', atomp_name
                        for k in range(2,len(finalT[i][j])):
                            exist1 = finalT[i][j].count(finalT[i][j][k])
                            exist2 = dicT[atomp_name].count(finalT[i][j][k])
                            if exist1 == 2*exist2:
                                continue
                            elif exist1 != 2*exist2:
                                expectexist = exist1/2
                                need_to_add = int(expectexist-exist2)
                                for l in range(need_to_add):
                                    dicT[atomp_name].append(finalT[i][j][k])
                    if len(dicT[atomp_name])!=((len(finalT[i][j])-2)/2):
                        print ('error!, did not paste enough data!', atomp_name)
    sum1 = 0
    for key in dicT:
        sum1 += len(dicT[key])
    
    for i in range(len(fNonbonding)):
        for j in range(len(fNonbonding[i])):
            for k in range(len(fNonbonding[i][j])):
                if fNonbonding[i][j][k][0][0] != 'H' and fNonbonding[i][j][k][2][0] != 'H':
                    atomp_name3 = fNonbonding[i][j][k][0]+'__'+fNonbonding[i][j][k][2]
                    atomp_name4 = fNonbonding[i][j][k][2]+'__'+fNonbonding[i][j][k][0]
                    ind1 = fNonbonding[i][j][k][1].index('_')
                    ind2 = fNonbonding[i][j][k][3].index('_')
                    r1 = fNonbonding[i][j][k][1][ind1+1:]
                    r2 = fNonbonding[i][j][k][3][ind2+1:]
                    ind3 = r1.index('_')
                    ind4 = r2.index('_')
                    resipair = r1[:ind3]+'_'+r2[:ind4]
                    if not dicN.has_key(atomp_name3) and not dicN.has_key(atomp_name4):
                        dicN[atomp_name3] = []
                        dicN[atomp_name3].append(fNonbonding[i][j][k][-1])
                        dicN_info[atomp_name3] = []
                        dicN_info[atomp_name3].append(resipair)
                    elif dicN.has_key(atomp_name3):
                        dicN[atomp_name3].append(fNonbonding[i][j][k][-1])
                        dicN_info[atomp_name3].append(resipair)
                    elif dicN.has_key(atomp_name4):
                        dicN[atomp_name4].append(fNonbonding[i][j][k][-1])
                        dicN_info[atomp_name4].append(resipair)
    #return fNonbonding
    return dicT, dicN, dicN_info


        

                        



        
