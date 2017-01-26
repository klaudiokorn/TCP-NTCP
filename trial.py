# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt

dic_general = {}
dic_structure = {}
D, d_DV, V_tot = [], [], []

########### RESCATE DE DATOS
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False

def n_structures(archivo):
    n=0
    for linea in archivo:
        if linea.startswith('Structure:'): 
            n+=1
    return n
    
def RescueDataGeneral(archivo):
    for linea in archivo:
        linea = linea.strip()
      
        part_general = linea.partition(':')
        dic_general[part_general[0].strip()]= part_general[2]
#    return dic_general

def RescueDataStructure(archivo):
    for linea in archivo:
        linea = linea.strip()
        
        part_structure = linea.partition(':')
        dic_structure[part_structure[0]]= part_structure[2]
#    return dic_structure

def RescueDataDoseVolume(archivo):
    q = []
    for linea in archivo:
        linea = linea.strip()
        q = linea.split()
        if len(q) == 3 and is_number(q[0]) == True: 
            D.append(float(q[0]))
            d_DV.append(float(q[2]))    



##### FUNCIONES PARA LAS DOSIS
def D_Reference(D_pres, curve_pres):
	D_ref = D_pres/(curve_pres/100.)
	return D_ref


def Norm_Dose(D, D_ref):
	D_norm = [0.]*len(D)
	for i in range(len(D)):
		D_norm[i] = D[i]/D_ref
	return D_norm	

def d_D(D):
	dD = (D[-1]-D[0])/(float(len(D)) - 1)
	return dD


##### FUNCIONES PARA LOS VOLUMENES

def Volume_total(d_DV, D, dD):
	V_tot = np.sum(d_DV)*dD
	return V_tot

def dif_Norm_Volume(d_DV, V_tot, dD):
	dV_norm = [0.]*len(d_DV)
	for i in range(len(d_DV)):
		dV_norm[i] = (d_DV[i]*dD)/V_tot
	return dV_norm	


def acc_Volume(d_DV, dD):
	acc_V = [0.]*len(d_DV)	
	for i in range(len(d_DV)):
		if i == 0: acc_V[i] = d_DV[0]*dD
		acc_V[i] = d_DV[i]*dD + acc_V[i-1]
	return acc_V	

def acc_Norm_Volume(dV_norm):
	acc_V_norm = [0.]*len(dV_norm)	
	for i in range(len(dV_norm)):
		if i == 0: acc_V_norm[i] = dV_norm[0]
		acc_V_norm[i] = dV_norm[i] + acc_V_norm[i-1]
	return acc_V_norm

def inv_acc_Norm_Volume(acc_V_norm):
	inv_acc_V_norm = [0.]*len(acc_V_norm)
	for i in range(len(acc_V_norm)):
		inv_acc_V_norm[i] = 1. - acc_V_norm[i]
	return inv_acc_V_norm


################### FUNCIONES DE PARAMETROS de DOSIS
n_fx = 25.      #fracciones
alfabeta = 15.  #alfa/beta
gamma = 4.      #gradiente dosis/respuestas
D_50 = 50.      #Para TCP (tumor)

SF2 = 0.5  #fraccion de supervivencia celular a Dref= 2 Gy
Dref = 2.  #Dosis ref EUD
m = 1      #efecto volumen 

def EUD_Niemierko(dV_norm, D, SF2, Dref):
	eud = [0.]*len(dV_norm)
	for i in range(len(dV_norm)):
		eud[i]= dV_norm[i]*SF2**(D[i]/Dref)

	EUD_niemierko= (Dref* np.log(np.sum(eud)))/(np.log(SF2))
	return EUD_niemierko


def EUD_Kutcher(d_DV, D, dD, m, V_tot, D_max):
	v = [0.]*len(d_DV)
	for i in range(len(d_DV)):
		v[i] = d_DV[i]*dD*(D[i]/D_max)**(1./m)
	V_eff = (np.sum(v)/V_tot)*100.
	EUD_kutcher = D_max*(V_eff/100.)**m
	return EUD_kutcher

############### TCP
def TCP_MeanDose(D_50, D_mean, gamma):
	TCP_D_mean = (1.+(D_50/D_mean)**(4.*gamma))**(-1.)	
	return TCP_D_mean

def TCP_EUD_Niemierko(D_50, EUD_niemierko, gamma):
	TCP_EUD_niemierko = (1.+(D_50/EUD_niemierko)**(4.*gamma))**(-1.)	
	return TCP_EUD_niemierko

	
################ NTCP

def NTCP_EUD_Kutcher(TD_50_V1, EUD_kutcher, gamma):
	NTCP_EUD_kutcher = (1.+(TD_50_V1/EUD_kutcher)**(4.*gamma))**(-1.)	
	return NTCP_EUD_kutcher


################# Equivalente 2 Gy/fx


def Equiv_2Gy_fx(alfabeta, n_fx, D):
	equiv_2Gy_fx = [0.]*len(D)
	for i in range(len(D)):
		equiv_2Gy_fx[i] = D[i]*(alfabeta + (D[i])/n_fx)/(alfabeta + 2.)
	return equiv_2Gy_fx

################# Dosis efectivas
#para ptv
def D_eff_Brahme(D_mean, D_STD, gamma, TCP_D_mean):
	D_eff_brahme = D_mean*(1. - (gamma/(2.*TCP_D_mean))*(D_STD/D_mean)**2 )
	return D_eff_brahme


#para OAR
def D_eff_Lyman(D, inv_acc_V_norm):
	D_eff_lyman = [0.]*len(D)
	D_rev = D[::-1]
	inv_acc_V_norm_rev = inv_acc_V_norm[::-1]	
#	print D_rev

	for i in range(len(D)):
		if i == 0:
			D_eff_lyman[0] = D_rev[0]
		else : D_eff_lyman[i] = D_rev[i] + (D_eff_lyman[i-1] - D_rev[i])*(inv_acc_V_norm_rev[i-1]/inv_acc_V_norm_rev[i]) 

#	print D_eff_lyman
	return D_eff_lyman[-1]

#Para hallar D_98%, D_50% y D_2%

def D_098(D, inv_acc_V_norm):
	inv_acc_V_norm_rev = inv_acc_V_norm[::-1]
	D_rev = D[::-1]		
	d_098 = np.interp(0.980, inv_acc_V_norm_rev, D_rev)
	return d_098

def D_002(D, inv_acc_V_norm):
	inv_acc_V_norm_rev = inv_acc_V_norm[::-1]
	D_rev = D[::-1]		
	d_002 = np.interp(0.020, inv_acc_V_norm_rev, D_rev)
	return d_002

def D_050(D, inv_acc_V_norm):
	inv_acc_V_norm_rev = inv_acc_V_norm[::-1]
	D_rev = D[::-1]		
	d_050 = np.interp(0.50, inv_acc_V_norm_rev, D_rev)
	return d_050

def HI(D, inv_acc_V_norm):
	hi = (D_002(D, inv_acc_V_norm) - D_098(D, inv_acc_V_norm))/D_050(D, inv_acc_V_norm)
	return hi


######### MAIN #########

name="paciente1PTV.txt"

with open(name, "r") as f1:
    n = n_structures(f1)

with open(name, "r") as f1:
    RescueDataGeneral(f1)

with open(name, "r") as f1:
    RescueDataStructure(f1)
    
with open(name, "r") as f1:
    RescueDataDoseVolume(f1)


dD    = d_D(D)
acc_V = acc_Volume(d_DV, dD)
#print dD     
#V_tot     = float(dic_structure['Volume [cm³]'])
D_pres     = float(dic_general['Prescribed dose [Gy]'])
curve_pres = float(dic_general['% for dose (%)'])
D_max      = float(dic_structure['Max Dose [Gy]'])
D_mean     = float(dic_structure['Mean Dose [Gy]'])
D_STD      = float(dic_structure['STD [Gy]'])
V_tot          = Volume_total(d_DV, D, dD)
dV_norm        = dif_Norm_Volume(d_DV, V_tot, dD)
D_ref          = D_Reference(D_pres, curve_pres)
D_norm         = Norm_Dose(D, D_ref)
acc_V_norm     = acc_Norm_Volume(dV_norm)
inv_acc_V_norm = inv_acc_Norm_Volume(acc_V_norm)

d_098 = D_098(D, inv_acc_V_norm)
d_050 = D_050(D, inv_acc_V_norm)
d_002 = D_002(D, inv_acc_V_norm)
hi    = HI(D, inv_acc_V_norm)

equiv_2Gy_fx      = Equiv_2Gy_fx(alfabeta, n_fx, D)
equiv_2Gy_fx_norm = Equiv_2Gy_fx(alfabeta, n_fx, D_norm)

EUD_niemierko = EUD_Niemierko(dV_norm, D, SF2, Dref)
EUD_kutcher   = EUD_Kutcher(d_DV, D, dD, m, V_tot, D_max)

TCP_D_mean        = TCP_MeanDose(D_50, D_mean, gamma)
TCP_EUD_niemierko = TCP_EUD_Niemierko(D_50, EUD_niemierko, gamma)
NTCP_EUD_kutcher  = NTCP_EUD_Kutcher(D_50, EUD_kutcher, gamma)


D_eff_brahme = D_eff_Brahme(D_mean, D_STD, gamma, TCP_D_mean)
D_eff_lyman  = D_eff_Lyman(D, inv_acc_V_norm)


print ('\n')
print ('DATOS GENERALES')
print ('Nombre: ', dic_general['﻿Patient Name'])
print ('ID: ', dic_general['Patient ID'])
print ('Tipo: ', dic_general['Type'])
print ('Estado del Plan: ', dic_general['Plan Status'])
print ('Fecha: ', dic_general['Date'])

print ('Dosis Prescrita [Gy]: ', D_pres)
print ('Curva Prescrita [%]: ', curve_pres)
print ('D_ref [Gy]: ', np.around(D_ref, decimals=3))


print ('\n')
print ('Datos especificos:')
print ('Estructura: ', dic_structure['Structure'])               

#print 'suma dV/dD [cc]: ', Volume_total(d_DV, D, dD)
#print 'Volumen [cc]: ', dic_structure['Volume [cm³]']               
#print 'Volumen acc_V [cc]: ', acc_V[-1]

print ('D_max  [Gy]:', np.around(D_max, decimals=3))
print ('D_mean [Gy]:', np.around(D_mean, decimals=3),' u$\pm$ ',np.around(D_STD, decimals=3))
#print 'D_STD [Gy]:', D_STD
print ('D_98%  [Gy]:', np.around(d_098, decimals=3))
print ('D_50%  [Gy]:', np.around(d_050, decimals=3))
print ('D_2%   [Gy]:', np.around(d_002, decimals=3))
print ('HI:', np.around(hi, decimals=4))

print ('\n')
print ('Aplica en caso Estructura = Tumor:')
print ('D_eff_Brahme  [Gy]: ', np.around(D_eff_brahme, decimals=3))
print ('EUD_Niemierko [Gy]: ', np.around(EUD_niemierko, decimals=3))
print ('TCP_D_mean        : ', np.around(TCP_D_mean, decimals=4))
print ('TCP_EUD_Niemierko : ', np.around(TCP_EUD_niemierko, decimals=4))

print ('\n')
print ('Aplica en caso Estructura = OAR:')
print ('D_eff_Lyman [Gy]: ', np.around(D_eff_lyman, decimals=3))
print ('EUD_Kutcher [Gy]: ', np.around(EUD_kutcher, decimals=3))
print ('NTCP_EUD_kutcher: ', np.around(NTCP_EUD_kutcher, decimals=4))
print ('\n')

#print (D)
########## PLOTEO ##########
fig = plt.figure()

ax1 = fig.add_subplot(221)
ax1.set_title(dic_structure['Structure'])
ax1.set_ylabel('dV/dD $[cm^{3}/Gy]$')
plt.plot (D, d_DV, 'b-', equiv_2Gy_fx, d_DV , 'g-')

ax2 = fig.add_subplot(222)
ax2.set_title(dic_structure['Structure'])
plt.plot (D_norm, d_DV, 'b-')

ax3 = fig.add_subplot(223)
ax3.set_xlabel('Dosis $[Gy]$')
ax3.set_ylabel('Vol')
plt.plot (D, inv_acc_V_norm, 'b-', D, dV_norm , 'r-')


ax4 = fig.add_subplot(224)
ax4.set_xlabel('Dosis')
#plt.plot (D_norm, dV_norm)
plt.plot (D_norm, inv_acc_V_norm, 'b-', D_norm, dV_norm , 'r-')


plt.show()

