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


# In[3]:


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


# In[4]:


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


# In[5]:


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


# In[6]:


vt=qload('vt')


# In[8]:


tlist=np.linspace(0,40,400)
trajn=500


# In[9]:


evolv1=[]
cn=200
Lindc=Lind(cn)
psi0=vt[0]
rho0=psi0*psi0.dag()
emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn)
opts = Options(rhs_reuse=True)
evolv1.append(emc)

for i in range(8):
    psi0=vt[i+1]
    rho0=psi0*psi0.dag()
    emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn, options=opts)
    evolv1.append(emc)


# In[10]:


qsave(evolv1,'evolv1')


# In[11]:


evolv2=[]
cn=400
Lindc=Lind(cn)
psi0=vt[0]
rho0=psi0*psi0.dag()
emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn)
opts = Options(rhs_reuse=True)
evolv2.append(emc)

for i in range(8):
    psi0=vt[i+1]
    rho0=psi0*psi0.dag()
    emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn, options=opts)
    evolv2.append(emc)


# In[12]:


qsave(evolv2,'evolv2')


# In[13]:


evolv3=[]
cn=600
Lindc=Lind(cn)
psi0=vt[0]
rho0=psi0*psi0.dag()
emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn)
opts = Options(rhs_reuse=True)
evolv3.append(emc)

for i in range(8):
    psi0=vt[i+1]
    rho0=psi0*psi0.dag()
    emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn, options=opts)
    evolv3.append(emc)


# In[14]:


qsave(evolv3,'evolv3')


# In[15]:


evolv4=[]
cn=800
Lindc=Lind(cn)
psi0=vt[0]
rho0=psi0*psi0.dag()
emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn)
opts = Options(rhs_reuse=True)
evolv4.append(emc)

for i in range(8):
    psi0=vt[i+1]
    rho0=psi0*psi0.dag()
    emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn, options=opts)
    evolv4.append(emc)


# In[16]:


qsave(evolv4,'evolv4')


# In[ ]:


evolv5=[]
cn=1000
Lindc=Lind(cn)
psi0=vt[0]
rho0=psi0*psi0.dag()
emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn)
opts = Options(rhs_reuse=True)
evolv5.append(emc)

for i in range(8):
    psi0=vt[i+1]
    rho0=psi0*psi0.dag()
    emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn, options=opts)
    evolv5.append(emc)


# In[ ]:


qsave(evolv5,'evolv5')


# In[11]:


evolv1=qload('evolv1')
evolv2=qload('evolv2')


# In[17]:


def f(x,A):
    return(np.exp(-A*x))


# In[18]:


drate=[]
emsf=[]
for i in range(9):
    As=op.curve_fit(f,tlist,evolv1[i].expect[0])[0]
    drate.append(As[0])
    emsf.append(f(tlist,As))


# In[19]:


evolveplot=evolv1


# In[24]:


x_major_locator = MultipleLocator(20)
y_major_locator = MultipleLocator(0.4)

fig, axes = plt.subplots(1,1)
#plt.yscale('log')
axes.plot(tlist, np.log(evolveplot[0].expect[0]),linewidth=2,color='blue',label= r'$scar states$')

axes.plot(tlist, np.log(evolveplot[1].expect[0]),linewidth=1.7,color='darkgray')
axes.plot(tlist, np.log(evolveplot[2].expect[0]),linewidth=1.7,color='darkgray')
axes.plot(tlist, np.log(evolveplot[3].expect[0]),linewidth=1.7,color='darkgray')
axes.plot(tlist, np.log(evolveplot[4].expect[0]),linewidth=1.7,color='darkgray')
axes.plot(tlist, np.log(evolveplot[5].expect[0]),linewidth=1.7,color='darkgray')
axes.plot(tlist, np.log(evolveplot[6].expect[0]),linewidth=1.7,color='darkgray')
axes.plot(tlist, np.log(evolveplot[7].expect[0]),linewidth=1.7,color='darkgray')
axes.plot(tlist, np.log(evolveplot[8].expect[0]),linewidth=1.7,color='darkgray',label= r'$thermal states$')
axes.plot(tlist, np.log(emsf[0]),linewidth=2,linestyle='-.',color='red',label= r'$fit e^{-\alpha_i t} for scar states$')
axes.plot(tlist, np.log(emsf[1]),linewidth=1.7,linestyle='-.',color='black')
axes.plot(tlist, np.log(emsf[2]),linewidth=1.7,linestyle='-.',color='black')
axes.plot(tlist, np.log(emsf[3]),linewidth=1.7,linestyle='-.',color='black')
axes.plot(tlist, np.log(emsf[4]),linewidth=1.7,linestyle='-.',color='black')
axes.plot(tlist, np.log(emsf[5]),linewidth=1.7,linestyle='-.',color='black')
axes.plot(tlist, np.log(emsf[6]),linewidth=1.7,linestyle='-.',color='black')
axes.plot(tlist, np.log(emsf[7]),linewidth=1.7,linestyle='-.',color='black')
axes.plot(tlist, np.log(emsf[8]),linewidth=1.7,linestyle='-.',color='black',label= r'$fit e^{-\alpha_i t} for thermal states$')

plt.xlim([0,40])
#plt.ylim([0.45,1])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
axes.xaxis.set_major_locator(x_major_locator)
axes.yaxis.set_major_locator(y_major_locator)
plt.xlabel('t',fontsize=20)
plt.ylabel(r'$|\langle E_i | \rho(t) |E_i \rangle|^2$',fontsize=20)
axes.legend(loc=1)
axes.legend(framealpha=0.5,fontsize=15)
plt.rc('axes',linewidth=2)
plt.savefig('spin1EigenDecay.svg')


# In[ ]:




