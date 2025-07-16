import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csr_matrix
from qutip import *
from matplotlib.ticker import MultipleLocator

N=9               #total spin number
j=1            #spin size=1, may be changed to spin-3/2,spin-2, etc
d=int(2*j+1)  

idd=qeye(d)         #single site identity

ids=idd             #entire identity
for i in range(N-1):
    ids=tensor(ids,idd)

Em=np.zeros((d,d))

em=Qobj(Em)         

ec=em               
for i in range(N-1):
    ec=tensor(ec,em)

Sx=np.zeros((d,d))
for i in range(d-1):
    m=j-i
    Sx[i,i+1]=((j+m)*(j-m+1))**(1/2)/2
    Sx[i+1,i]=((j+m)*(j-m+1))**(1/2)/2
sx=Qobj(Sx)         #s_x

Su=np.zeros((d,d)) 
Su[1,0]=(j/2)**(1/2)
Su[0,1]=(j/2)**(1/2)
su=Qobj(Su)

Sd=np.zeros((d,d))  
Sd[d-2,d-1]=(j/2)**(1/2)
Sd[d-1,d-2]=(j/2)**(1/2)
sd=Qobj(Sd)

Lu=np.zeros((d,d)) 
Lu[0,1]=1
lu=Qobj(Lu)

Ld=np.zeros((d,d))  
Ld[d-1,d-2]=1
ld=Qobj(Ld)

Pu=np.zeros((d,d)) 
Pu[0,0]=1
pu=Qobj(Pu)

Pd=np.zeros((d,d))  
Pd[d-1,d-1]=1
pd=Qobj(Pd)

H0=ec               #H_0
for i in range(N):
    if i==0:
        h0=sx
    else:
        h0=idd
    for k in range(N-1):
        if k+1==i:
            h0=tensor(h0,sx)
        else:
            h0=tensor(h0,idd)
    H0=H0+h0

H=H0               #Hamiltonian for scar model
for i in range(N):
    if i==0:
        h1=pu
        h2=su
    elif i==N-1:
        h1=sd
        h2=pd
    else:
        h1=idd
        h2=idd
    for k in range(N-1):
        if k+1==i:
            h1=tensor(h1,pu)
            h2=tensor(h2,su)
        elif k==i:
            h1=tensor(h1,sd)
            h2=tensor(h2,pd)
        else:
            h1=tensor(h1,idd)
            h2=tensor(h2,idd)
    H=H-h1-h2
    
Pb=ids
for i in range(N):
    if i==0:
        pp=pu
    elif i==N-1:
        pp=pd
    else:
        pp=idd
    for k in range(N-1):
        if k+1==i:
            pp=tensor(pp,pu)
        elif k==i:
            pp=tensor(pp,pd)
        else:
            pp=tensor(pp,idd)
    Pb=Pb*(ids-pp)
    
Hp=Pb*H0*Pb


def Lind(c):
    gamma=((2*j)**(1/2))*c
    L=[]                #lindbladian operators for L
    for i in range(N):
        if i==0:
            sigma1=pu
            sigma2=lu
            sigma3=pu
        elif i==N-1:
            sigma1=pd
            sigma2=pd
            sigma3=ld
        else:
            sigma1=idd
            sigma2=idd
            sigma3=idd
        for k in range(N-1):
            if k+1==i:
                sigma1=tensor(sigma1,pu)
                sigma2=tensor(sigma2,lu)
                sigma3=tensor(sigma3,pu)
            elif k==i:
                sigma1=tensor(sigma1,pd)
                sigma2=tensor(sigma2,pd)
                sigma3=tensor(sigma3,ld)
            else:
                sigma1=tensor(sigma1,idd)
                sigma2=tensor(sigma2,idd)
                sigma3=tensor(sigma3,idd)
        ldb=(gamma**(1/2))*(sigma1-(1j/c)*sigma2-(1j/c)*sigma3)
        L.append(ldb)
    return(L)

def Linda(c):
    gamma=((2*j)**(1/2))/c
    L=[]                #lindblad operators for L'
    for i in range(N):
        if i==0:
            sigma2=lu
            sigma3=pu
        elif i==N-1:
            sigma2=pd
            sigma3=ld
        else:
            sigma2=idd
            sigma3=idd
        for k in range(N-1):
            if k+1==i:
                sigma2=tensor(sigma2,lu)
                sigma3=tensor(sigma3,pu)
            elif k==i:
                sigma2=tensor(sigma2,pd)
                sigma3=tensor(sigma3,ld)
            else:
                sigma2=tensor(sigma2,idd)
                sigma3=tensor(sigma3,idd)
        ldb=(gamma**(1/2))*(sigma2+sigma3)
        L.append(ldb)
    return(L)



st1=basis(d,0)
st2=basis(d,1)
st3=basis(d,2)

psi0=st1           #initial states
for i in range(N-1):
    psi0=tensor(psi0,st1)
    
psi1=tensor(st1,st1,st1,st2,st2,st1,st1,st1,st1)
psi2=tensor(st1,st1,st2,st2,st3,st2,st3,st2,st1)


tlist=np.linspace(0,30,300)

#calculation for evolution
trajn=1000

oper=psi0*psi0.dag()
evsc=mesolve(H, psi0, tlist,[], [])
evL1=mcsolve(H0, psi0, tlist,Lind(50), [], ntraj=trajn)
evL2=mcsolve(H0, psi0, tlist,Lind(100), [], ntraj=trajn)
evL3=mcsolve(H0, psi0, tlist,Lind(200), [], ntraj=trajn)

evLa1=mcsolve(Hp, psi0, tlist,Linda(50), [], ntraj=trajn)
evLa2=mcsolve(Hp, psi0, tlist,Linda(100), [], ntraj=trajn)
evLa3=mcsolve(Hp, psi0, tlist,Linda(200), [], ntraj=trajn)

#calculation for evolution under other initial states

oper2=psi2*psi2.dag()
evsc1=mesolve(H, psi2, tlist,[], [])
evL11=mcsolve(H0, psi2, tlist,Lind(50), [], ntraj=trajn)
evL21=mcsolve(H0, psi2, tlist,Lind(100), [], ntraj=trajn)
evL31=mcsolve(H0, psi2, tlist,Lind(200), [], ntraj=trajn)

evLa11=mcsolve(Hp, psi2, tlist,Linda(50), [], ntraj=trajn)
evLa21=mcsolve(Hp, psi2, tlist,Linda(100), [], ntraj=trajn)
evLa31=mcsolve(Hp, psi2, tlist,Linda(200), [], ntraj=trajn)

#overlap with initial states
evscE=expect(psi0*psi0.dag(),evsc.states)
evL1E=expect(psi0*psi0.dag(),evL1.states)
evL2E=expect(psi0*psi0.dag(),evL2.states)
evL3E=expect(psi0*psi0.dag(),evL3.states)
evLa1E=expect(psi0*psi0.dag(),evLa1.states)
evLa2E=expect(psi0*psi0.dag(),evLa2.states)
evLa3E=expect(psi0*psi0.dag(),evLa3.states)

x_major_locator = MultipleLocator(10)
y_major_locator = MultipleLocator(0.5)

fig, axes = plt.subplots(1,1)

axes.plot(tlist, evscE,linewidth=2, label=r'$H$',color='b')
axes.plot(tlist, evL1E,linewidth=2, label=r'$\mathcal{L}_{+},c=50$',color='silver')
axes.plot(tlist, evL2E,linewidth=2, label=r'$\mathcal{L}_{+},c=100$',color='dimgrey')
axes.plot(tlist, evL3E,linewidth=2, label=r'$\mathcal{L}_{+},c=200$',color='black')
axes.plot(tlist, evLa1E,linewidth=2, linestyle='-.', label=r'$\mathcal{L}^{’}_{+},c=50$',color='darkred')
axes.plot(tlist, evLa2E,linewidth=2, linestyle='-.', label=r'$\mathcal{L}^{’}_{+},c=100$',color='red')
axes.plot(tlist, evLa3E,linewidth=2, linestyle='-.', label=r'$\mathcal{L}^{’}_{+},c=200$',color='salmon')
#plt.xlim(12.55,12.75)
#plt.ylim(0.76,0.79)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
axes.xaxis.set_major_locator(x_major_locator)
axes.yaxis.set_major_locator(y_major_locator)
plt.xlabel('t',fontsize=20)
plt.ylabel(r'$|\langle \psi(0) | \rho(t) |\psi(0) \rangle|^2$',fontsize=20)
#axes.legend(loc=1)
#axes.legend(framealpha=0.5)
plt.rc('axes',linewidth=2)
plt.savefig('spin1psi0.svg')

x_major_locator = MultipleLocator(0.5)
y_major_locator = MultipleLocator(0.2)

fig, axes = plt.subplots(1,1)

axes.plot(tlist, evscE,linewidth=2, label=r'$H$',color='b')
axes.plot(tlist, evL1E,linewidth=2, label=r'$\mathcal{L}_{+},c=50$',color='silver')
axes.plot(tlist, evL2E,linewidth=2, label=r'$\mathcal{L}_{+},c=100$',color='dimgrey')
axes.plot(tlist, evL3E,linewidth=2, label=r'$\mathcal{L}_{+},c=200$',color='black')
axes.plot(tlist, evLa1E,linewidth=2, linestyle='-.', label=r'$\mathcal{L}^{’}_{+},c=50$',color='darkred')
axes.plot(tlist, evLa2E,linewidth=2, linestyle='-.', label=r'$\mathcal{L}^{’}_{+},c=100$',color='red')
axes.plot(tlist, evLa3E,linewidth=2, linestyle='-.', label=r'$\mathcal{L}^{’}_{+},c=200$',color='salmon')
plt.xlim(19.5,20.6)
plt.ylim(0.3,0.8)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
axes.xaxis.set_major_locator(x_major_locator)
axes.yaxis.set_major_locator(y_major_locator)
plt.xlabel('t',fontsize=20)
plt.ylabel(r'$|\langle \psi(0) | \rho(t) |\psi(0) \rangle|^2$',fontsize=20)
axes.legend(framealpha=0.5)
plt.rc('axes',linewidth=2)
axes.legend(prop={'size':15})
plt.savefig('spin1psi0sub.svg')

