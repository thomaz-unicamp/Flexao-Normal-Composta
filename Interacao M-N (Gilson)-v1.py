# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 14:09:23 2023

@author: Prof. Thomaz Buttignol
         Universidade Estadual de Campinas (Unicamp)
         Departamento de Estruturas (DES)
         Faculdade de Engenharia Civil, Arquitetura e Urbanismo (FECFAU)

Flexão Normal Composta
Interação M-N
Método do prof. Gilson B. Fernandes, DES-FECFAU, Unicamp
"""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

########### 1. Dados de entrada ###########

#forças atuantes:
print("")
print("Forças atuantes:")
Mk = float(input("Momento característico [kN.cm]:"))
Nk = float(input("Força característica [kN]:"))
gf = 1.4 # majoração dos esforços

#seção transversal 
print("")
print("seção transversal:")
bb = float(input("largura da seção [cm]:"))
hh = float(input("altura da seção [cm]:"))
dd = float(input("altura útil [cm]:"))

#materiais
print("")
print("materiais:")
fck = float(input("resistência à compressão [MPa]:"))
ecu = float(input("deformação última do concreto (ecu) [-]:"))
ec2 = float(input("deformação do concreto (ec2) [-]:"))
lmb = float(input("lambda [-]:"))
ac = float(input("alpha,c [-]:"))

fyk = 500 # em MPa - resistência característica do aço CA-50
Es = 210000 # em MPa - módulo de elasticidade do aço
esu = 0.01

########### 2. Inicialização ############

#parâmetros geométricos
dl = hh - dd
eta = dl/dd

#propriedades mecânicas
fcd = fck/1.4/10 # em kN/cm2 - resistência à compressão de cálculo
k = (ecu - ec2)/ecu
fyd = fyk/1.15 # em kN/cm2 - resistência de cálculo do aço CA-50 

#forças externas
Nd = gf*Nk # em kN - força normal de cálculo

#excentricidade
e1 = Mk/Nk # em cm - excentricidade em relação ao CG
es = e1 + (dd - dl)/2 # em cm - excenrtricidade em relação a As

#esforços adimensionais
msd = Nd*es/(bb*dd**2*fcd) # momento adimensional de cálculo
nsd = Nd/(bb*dd*fcd) # força cortante adimensional de cálculo

#linha neutra
bx11 = 1/lmb-1/lmb*np.sqrt(1-2*msd/ac)
bx12 = 1/lmb+1/lmb*np.sqrt(1-2*msd/ac)
bx21 = eta/lmb-eta/lmb*np.sqrt(1+2/(eta**2*ac)*((1-eta)*nsd-msd))
bx22 = eta/lmb+eta/lmb*np.sqrt(1+2/(eta**2*ac)*((1-eta)*nsd-msd))

########### 3. Solução do gráfico mi x bx ###########

#pontos de interesse
bA = eta
bB = 1
bC = 1/lmb
bD = (1+eta)/lmb
bE = 2
miA = lmb*ac*eta*(1-lmb*eta/2)
miB = lmb*ac*(1-lmb/2)
miC = ac/2
miD = ac/2*(1-eta**2)
miE = miD

pt_mi = pd.DataFrame({"beta":[0, bA, bB, bC, bD, bE], "mi":[0, miA, miB, miC, miD, miE]}) #dataframe

#gráfico mi x bx
beta = []
mi = []
steps = 100
for i in range (steps+1):    
    if i == 0: 
        beta.append(0)
        mi.append(0)
    if i < steps:
        beta.append(bD/steps*(i+1))
        mi.append(lmb*ac*beta[i]*(1-lmb/2*beta[i]))
    else: 
        beta.append(2)
        mi.append(ac/2*(1-eta**2))

values = {"beta":beta, "mi":mi}
gf_mi = pd.DataFrame(values)
fig, ax = plt.subplots(ncols=2)
sns.lineplot(data=values, x="beta", y="mi", color = 'black', ax = ax[0])

#plot
pts = pt_mi.to_numpy()
pts = np.transpose(pts)
sns.lineplot(x=[0,bE], y=[msd,msd],  color = 'red', ax = ax[0])

#cálculo de bx1
if msd > miC: bx1 = 'NaN'
else:
    for i in range (len(pts[1])):
        if pts[1,i-1] < msd < pts[1,i]:
            if pts[0,i-1] < bx11 < pts[0,i]:
                bx1 = bx11
            else: bx1 = bx12        
print("")
print("bx1 =", bx1)
#intervalo de msd (ver artigo "Dimensionamento de seções de concreto armado submetidos a flexão composta normal")
if 0 < msd <= miA: inf1 = round(bx1,3); sup1 = round(bA,3); print("intervalo 1 =>", inf1, "<= bx <", sup1) #Caso A
if miA < msd <= miB: inf1 = round(bA,3); sup1 = round(bx1,3); print("intervalo 1 =>", inf1, "< bx <=", sup1) #Caso B
if miB < msd <= miC: inf1 = 1 ; sup1 = round(bx1,3);print("intervalo 1 =>", inf1, "< bx <=", sup1) #Caso C1
if miC < msd <= miD: inf1 = round(bD,3); sup1 = round(bx1,3); print("intervalo 1 =>", inf1, "< bx <=", sup1) #Caso C2
if msd > miD: inf1 = round(bA,3); sup1 = 1E10; print("intervalo 1 =>","bx >", inf1) #Caso D (não existe bx1)

########### 4. Solução do gráfico ni x bx ###########

#pontos de interesse
bF = eta
bG = eta/lmb
bH = 1
bI = (1+eta)/lmb
bJ = 2
niF = (2*msd-ac*lmb*eta**2*(2-lmb))/(2*(1-eta))
niG = (2*msd-ac*eta**2)/(2*(1-eta))
niH = (2*msd-ac*lmb*(2*eta-lmb))/(2*(1-eta))
niI = (2*msd+ac*(1-eta**2))/(2*(1-eta))
niJ = niI

pt_ni = pd.DataFrame({"beta":[0, bF, bG, bH, bI, bJ], "ni":[msd/(1-eta), niF, niG, niH, niI, niJ]}) #dataframe

#gráfico ni x bx
beta = []
ni = []
steps = 100
for i in range (steps+1):    
    if i == 0: 
        beta.append(0)
        ni.append(msd/(1-eta))
    if i < steps:
        beta.append(bI/steps*(i+1))
        ni.append(lmb*ac*beta[i]+(msd-lmb*ac*beta[i]*(1-lmb/2*beta[i]))/(1-eta))
    else: 
        beta.append(2)
        ni.append(ac*(1+eta)+(msd-ac/2*(1-eta**2))/(1-eta))

values = {"beta":beta, "ni":ni}
gf_ni = pd.DataFrame(values)
sns.lineplot(data=values, x="beta", y="ni", color = 'black', ax = ax[1])

#plot
pts = pt_ni.to_numpy()
pts = np.transpose(pts)
plt.plot([0,bJ],[nsd,nsd], color = 'red')
plt.show()

#cálculo de bx2
if (nsd < niG) or (nsd > niI): bx2 = 'NaN'
else:
    for i in range (len(pts[1])):
        if pts[1,i-1] < nsd < pts[1,i]:
            if pts[0,i-1] < bx21 < pts[0,i]:
                bx2 = bx21
            else: bx2 = bx22        
print("")
print("bx2 =", bx2)

#intervalo de nsd (ver artigo "Dimensionamento de seções de concreto armado submetidos a flexão composta normal")
if nsd < niG: inf2 = 0; sup2 = 1E10; print("intervalo 2 =>",inf2, "< bx <", sup2) #Caso A (não existe bx2)
if niG < nsd <= niI: 
    if bx2 < 1: inf2 = round(bx2,3); sup2 = 1; print("intervalo 2 =>", inf2, "<= bx <", sup2); #Caso 2A
    if bx2 == 1: inf2 = 1; sup2 = 1; print("intervalo 2 =>", inf2, "= bx =", sup2) #Caso 2B
    if bx2 > 1: inf2 = 1; sup2 = round(bx2,3); print("intervalo 2 =>", inf2, "< bx <=", sup2) #Caso 2C
if nsd > niI: inf2 = 1; sup2 = 1E10; print("intervalo 2 =>",inf2, "< bx <", sup2) #Caso A (não existe bx2)

##### Solução #####
print("")
inf = max (inf1, inf2)
sup = min(sup1, sup2)

if sup1 < inf2:
    print("solução =>","bx2", ">", "bx1")
else:
    print("solução =>", inf, "< bx <", sup)


########## 5. Cálculo da armadura ############

#entrada
print("")
bx = float(input("Entre com o valor de bx:")) #valor de entrada
print(f'bx de cálculo = {bx}')

#parâmetros adimensionais (omega e mi)
if bD > bx:
    w = lmb*ac*bx #omega (força resistente do concreto adimensionalizada)
    mi = lmb*ac*bx*(1-lmb/2*bx)
else:
    w = ac*(1+eta)
    mi = ac/2*(1-eta**2)

#domínios 2, 3, 4 e 4a
if bx == eta: Ass = 0
if (0 < bx < eta) or (eta < bx <= 1 + eta):
    sigss = Es/10*ecu/bx*(bx-eta)
    if sigss > 43.5: sigss = 43.5   

#domínio 5
if bx > 1 + eta:
    ec2/(bx - k*(1+eta))

if bx < bA: Ass = -(msd - mi)/(1-eta)*bb*dd*fcd/sigss
else: Ass = (msd - mi)/(1-eta)*bb*dd*fcd/sigss

print("Área de aço:")
print("Ass=", round(Ass,2), "cm2")

#cálculo de As
#domínio 2a
if bx <= eta: 
    sigs = 43.5
    As = -(Nd - bb*dd*lmb*bx*ac*fcd + Ass*sigss)/sigs

#domínios 2b e 3
if eta  < bx <= 0.63: 
    sigs = 43.5
    As = -(Nd - bb*dd*lmb*bx*ac*fcd - Ass*sigss)/sigs

#domínio 4
if 0.63 < bx < 1:
    sigs = Es/10*ecu/bx*(1-bx)
    As = -(Nd - bb*dd*lmb*bx*ac*fcd - Ass*sigss)/sigs

if bx == 1:
    sigs = 0
    As = 0

#domínio 4a
if 1 < bx <= 1+eta:
    sigs = Es/10*ecu/bx*(bx-1)
    As = (Nd - bb*dd*lmb*bx*ac*fcd - Ass*sigss)/sigs

#domínio 5 (caso 1)
if 1+eta < bx <= (1+eta)/lmb:
    sigs = Es/10*ec2/(bx-k*(1+eta))*(bx-1)
    As = (Nd - bb*dd*lmb*bx*ac*fcd - Ass*sigss)/sigs

#domínio 5 (caso 2)
if bx > (1+eta)/lmb:
    sigs = Es/10*ec2/(bx-k*(1+eta))*(bx-1)
    As = (Nd - bb*hh*ac*fcd - Ass*sigss)/sigs

print("As=", round(As,2), "cm2")

########## 6. Solução para par (bx, A) ##########

print("")
var = float(input("Entre com o valor da relação As/Ass:")) #valor de entrada
ns = var/2
nss = 1 - ns
print("As =",ns,"A","e Ass =", nss,"A")

from scipy.optimize import fsolve

def equations(vars):
    bx, A = vars
    
    #tensões nas armaduras para os domínios de deformação
    if bx < eta: #domínio 2a
        sigs = 43,5; sigss = Es/10*ecu/bx*(bx-eta) 
        eq1 = Nd - bb*dd*lmb*bx*ac*fcd + A*nss*sigss + A*ns*sigs
        eq2 = Nd*es - bb*dd**2*lmb*bx*ac*fcd*(1-0.5*lmb*bx) + A*nss*sigss*(dd-dl)        
    
    if bx == eta: 
        sigss = 0; sigs = 43.5
        eq1 = Nd - bb*dd*lmb*bx*ac*fcd + A*ns*sigs
        eq2 = Nd*es - bb*dd**2*lmb*bx*ac*fcd*(1-0.5*lmb*bx)   
    
    if eta < bx < 1: #domínios 2, 3 e 4
        sigs = Es/10*ecu/bx*(1-bx)
        sigss = Es/10*ecu/bx*(bx-eta)
        if sigss > 43.5: sigss = 43.5
        if sigs > 43.5: sigs = 43.5
        eq1 = Nd - bb*dd*lmb*bx*ac*fcd - A*nss*sigss + A*ns*sigs
        eq2 = Nd*es - bb*dd**2*lmb*bx*ac*fcd*(1-0.5*lmb*bx) - A*nss*sigss*(dd-dl)      
    
    if bx == 1: 
        sigs = 0; sigss = 43.5
        eq1 = Nd - bb*dd*lmb*bx*ac*fcd - A*nss*sigss
        eq2 = Nd*es - bb*dd**2*lmb*bx*ac*fcd*(1-0.5*lmb*bx) - A*nss*sigss*(dd-dl) 
    
    if 1 < bx <= 1 + eta: #domínio 4a
        sigs = Es/10*ecu/bx*(bx-1) ; sigss = 43.5       
        eq1 = Nd - bb*dd*lmb*bx*ac*fcd - A*nss*sigss - A*ns*sigs
        eq2 = Nd*es - bb*dd**2*lmb*bx*ac*fcd*(1-0.5*lmb*bx) - A*nss*sigss*(dd-dl) 
        
    if bx > 1 + eta: #domínio 5
        sigs = Es/10*ec2/(bx-k*(1+eta))*(bx-1); sigss = 43.5
        eq1 = Nd - bb*dd*lmb*bx*ac*fcd - A*nss*sigss - A*ns*sigs
        eq2 = Nd*es - bb*dd**2*lmb*bx*ac*fcd*(1-0.5*lmb*bx) - A*nss*sigss*(dd-dl) 

    return [eq1, eq2]

ini = (1, 1)
bx, A =  fsolve(equations, ini)

print("")
print("bx =",round(bx,5))
print("As =",round(ns*A,2),"cm2", "e Ass =",round(nss*A,2),"cm2")
