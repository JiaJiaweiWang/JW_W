#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from scipy import optimize as op
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


# In[2]:


N=7               #链长
j=1            #自旋大小
d=int(2*j+1)      #自旋矩阵维度


# In[3]:


idd=qeye(d)         #单自旋单位算符

ids=idd             #集体自旋单位算符
for i in range(N-1):
    ids=tensor(ids,idd)

Em=np.zeros((d,d))

em=Qobj(Em)         #单格点空算符

ec=em               #整链空算符
for i in range(N-1):
    ec=tensor(ec,em)

Sx=np.zeros((d,d))
for i in range(d-1):
    m=j-i
    Sx[i,i+1]=((j+m)*(j-m+1))**(1/2)/2
    Sx[i+1,i]=((j+m)*(j-m+1))**(1/2)/2
sx=Qobj(Sx)         #单自旋

Su=np.zeros((d,d))  #最高自旋跃迁算符
Su[1,0]=(j/2)**(1/2)
Su[0,1]=(j/2)**(1/2)
su=Qobj(Su)

Sd=np.zeros((d,d))  #最低自旋跃迁算符
Sd[d-2,d-1]=(j/2)**(1/2)
Sd[d-1,d-2]=(j/2)**(1/2)
sd=Qobj(Sd)

Lu=np.zeros((d,d))  #最高耗散跃迁算符
Lu[0,1]=1
lu=Qobj(Lu)

Ld=np.zeros((d,d))  #最低耗散跃迁算符
Ld[d-1,d-2]=1
ld=Qobj(Ld)

Pu=np.zeros((d,d))  #最高自旋投影
Pu[0,0]=1
pu=Qobj(Pu)

Pd=np.zeros((d,d))  #最低自旋投影
Pd[d-1,d-1]=1
pd=Qobj(Pd)


# In[4]:


H0=ec               #哈密顿量
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


# In[5]:


def Lind(c):
    gamma=((2*j)**(1/2))*c
    L=[]                #lindblad算符
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

def Lind1(c):
    gamma=((2*j)**(1/2))/c
    L=[]                #lindblad算符
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


# In[6]:


H=H0               #疤痕哈密顿量
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


# In[7]:


def st(i):
    state=basis(d,i)
    return(state)

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
Ut=ec
for i1 in range(d):
    for i2 in range(d):
        for i3 in range(d):
            for i4 in range(d):
                for i5 in range(d):
                    for i6 in range(d):
                        for i7 in range(d):
                            Ut=Ut+tensor(st(i2),st(i3),st(i4),st(i5),st(i6),st(i7),st(i1))*tensor(st(i1).dag(),st(i2).dag(),st(i3).dag(),st(i4).dag(),st(i5).dag(),st(i6).dag(),st(i7).dag())

Pt=(ids+Ut+(Ut**2)+(Ut**3)+(Ut**4)+(Ut**5)+(Ut**6))/N


# In[8]:


Hp=Pb*Pt*H*Pt*Pb


# In[9]:


tlist=np.linspace(0,60,400)
trajn=200


# In[10]:


py=(tensor(st(0),st(1))+tensor(st(1),st(2)))/(2**(1/2))
Py=py*py.dag()
for i in range(N-2):
    Py=tensor(Py,idd)
Ny=Py
Pym=Py
for i in range(N-1):
    Pym=Ut.dag()*Pym*Ut
    Ny=Ny+Pym


# In[11]:


psi0=[]
phi1=tensor(st(0),st(0),st(0),st(0),st(0),st(0),st(0))
phi2=tensor(st(0),st(1),st(0),st(1),st(0),st(0),st(0))
phi3=tensor(st(0),st(1),st(0),st(1),st(0),st(1),st(0))
phi4=tensor(st(0),st(1),st(2),st(1),st(0),st(1),st(0))
phi5=tensor(st(1),st(1),st(1),st(1),st(1),st(1),st(1))

psi0.append(phi1)
psi0.append(phi2)
psi0.append(phi3)
psi0.append(phi4)
psi0.append(phi5)
psi0[0].dag()*Pb*psi0[0]


# In[ ]:


evolve=[]
cn=1200
Lindc=Lind(cn)
emc=mcsolve(H0, psi0[0], tlist,Lindc, [Pb], ntraj=trajn)
opts = Options(rhs_reuse=True)
evolve.append(emc)

for i in range(4):
    emc=mcsolve(H0, psi0[i+1], tlist,Lindc, [Pb], ntraj=trajn, options=opts)
    evolve.append(emc)


# In[244]:


fig, axes = plt.subplots(1,1)

axes.plot(tlist, evolve[0].expect[0],linewidth=2,color='green',label= r'$scar states$')
axes.plot(tlist, evolve[1].expect[0],linewidth=2,color='yellow')
axes.plot(tlist, evolve[2].expect[0],linewidth=2,color='silver',label= r'$scar states$')
axes.plot(tlist, evolve[3].expect[0],linewidth=1.5,color='deepskyblue')
axes.plot(tlist, evolve[4].expect[0],linewidth=1.5,color='deepskyblue')
#plt.xlim([0,40])
#plt.ylim([0.4,1])
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#axes.xaxis.set_major_locator(x_major_locator)
#axes.yaxis.set_major_locator(y_major_locator)
#plt.xlabel('t',fontsize=20)
#plt.ylabel(r'$|\langle E_i | \rho(t) |E_i \rangle|^2$',fontsize=20)
#axes.legend(loc=1)
#axes.legend(framealpha=0.5,fontsize=15)
#plt.rc('axes',linewidth=2)
#plt.savefig('spin1EigenDecay3.svg')


# In[245]:


Nyexp=[0.434349483238014,0.966544883597522,1.213875098100984,1.288186996849070,1.376855937373432]


# In[269]:


x_major_locator = MultipleLocator(20)
y_major_locator = MultipleLocator(0.1)

fig, axes = plt.subplots(1,1)

axes.plot(tlist, evolve[0].expect[0],linewidth=2.5,color='blue',label= r'$tr(\rho(t)\hat{P})$')
axes.plot(tlist, evolve[1].expect[0],linewidth=2.5,color='silver',label= r'$thermal$')
axes.plot(tlist, evolve[2].expect[0],linewidth=2.5,color='silver')
#axes.plot(tlist, evolve[3].expect[0],linewidth=1.5,color='green')
#axes.plot(tlist, evolve[4].expect[0],linewidth=1.5,color='blue')
axes.plot(tlist, np.exp(-2*((2*j)**(1/2))*Nyexp[0]/cn*tlist) ,linestyle='-.',linewidth=2.5,color='red',label= r'$e^{\gamma_2 \bar{N}_i t}$')
axes.plot(tlist, np.exp(-2*((2*j)**(1/2))*Nyexp[1]/cn*tlist) ,linestyle='-.',linewidth=2.5,color='black',label= r'$e^{\gamma_2 \bar{N}_i t}$')
axes.plot(tlist, np.exp(-2*((2*j)**(1/2))*Nyexp[2]/cn*tlist) ,linestyle='-.',linewidth=2.5,color='black')
#axes.plot(tlist, np.exp(-2*((2*j)**(1/2))*Nyexp[3]/cn*tlist) ,linestyle='-.',linewidth=1.5,color='silver')
#axes.plot(tlist, np.exp(-2*((2*j)**(1/2))*Nyexp[4]/cn*tlist) ,linestyle='-.',linewidth=1.5,color='silver')
plt.xlim([0,60])
plt.ylim([0.82,1])
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
axes.xaxis.set_major_locator(x_major_locator)
axes.yaxis.set_major_locator(y_major_locator)
plt.xlabel('t',fontsize=20)
plt.ylabel(r'$tr(\rho(t)\hat{P})$',fontsize=20)
axes.legend(loc=1)
axes.legend(framealpha=0.5,fontsize=15)
plt.rc('axes',linewidth=2)
plt.savefig('spin1InitialDecay.svg')


# In[222]:


def f(x,A):
    return(np.exp(-A*x))

drate=[]
emsf=[]
for i in range(5):
    As=op.curve_fit(f,tlist,evolve[i].expect[0])[0]
    drate.append(As[0])
    emsf.append(f(tlist,As))


# In[223]:


x_major_locator = MultipleLocator(100)
y_major_locator = MultipleLocator(0.2)

fig, axes = plt.subplots(1,1)

axes.plot(tlist, evolve[0].expect[0],linewidth=2,color='green',label= r'$tr(\rho(t)\hat{P})$')
axes.plot(tlist, evolve[1].expect[0],linewidth=2,color='green')
axes.plot(tlist, evolve[2].expect[0],linewidth=2,color='green')
#axes.plot(tlist, evolve[3].expect[0],linewidth=1.5,color='green')
#axes.plot(tlist, evolve[4].expect[0],linewidth=1.5,color='blue')
axes.plot(tlist, emsf[0] ,linestyle='-.',linewidth=2,color='silver',label= r'$e^{\gamma_2 \bar{N}_i t}$')
axes.plot(tlist, emsf[1] ,linestyle='-.',linewidth=2,color='silver')
axes.plot(tlist, emsf[2] ,linestyle='-.',linewidth=2,color='silver')
#axes.plot(tlist, emsf[3] ,linestyle='-.',linewidth=1.5,color='silver')
#axes.plot(tlist, emsf[4] ,linestyle='-.',linewidth=1.5,color='silver')
plt.xlim([0,200])
plt.ylim([0.5,1])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
axes.xaxis.set_major_locator(x_major_locator)
axes.yaxis.set_major_locator(y_major_locator)
plt.xlabel('t',fontsize=20)
plt.ylabel(r'$tr(\rho(t)\hat{P})$',fontsize=20)
axes.legend(loc=1)
axes.legend(framealpha=0.5,fontsize=15)
plt.rc('axes',linewidth=2)
plt.savefig('spin1EigenDecay3.png')





