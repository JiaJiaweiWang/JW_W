#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from matplotlib.ticker import MultipleLocator


# In[2]:


N=9               #链长
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


# In[5]:


st1=basis(d,0)
st2=basis(d,1)
st3=basis(d,2)
sz=Qobj([[1,0,0],[0,0,0],[0,0,-1]])
sx=Qobj([[0,1,0],[1,0,1],[0,1,0]])
psi0=tensor(st1,st1,st1,st1,st1,st1,st1,st1,st1)           #初态
psi1=tensor(st1,st2,st3,st2,st1,st1,st1,st1,st1)           #初态

rdd1=idd
for i in range(N-2):
    rdd1=tensor(rdd1,idd)

rdd2=idd
for i in range(N-3):
    rdd2=tensor(rdd2,idd)

s2=(tensor(st1,st1)-tensor(st3,st3))/(2**(1/2))
S2z=tensor(s2*s2.dag(),rdd2)


# In[6]:


tlist=np.linspace(0,120,600)


# In[7]:


evol0=mesolve(H, psi0, tlist,[], [psi0*psi0.dag(),S2z])


# In[8]:


fig, axes = plt.subplots(1,1)

axes.plot(tlist, evol0.expect[0],linewidth=2.0, label=r'$initialstates$',color='silver')
axes.plot(tlist, evol0.expect[1],linewidth=1.0, label=r'$obser2$')
#plt.xlim(19,21)
#plt.ylim(0.3,0.8)
axes.legend(loc=1);


# In[127]:


Ocan=np.zeros(len(tlist))
for i in range(len(tlist)):
    Ocan[i]=0.106169698658641
Ogra=np.zeros(len(tlist))
for i in range(len(tlist)):
    Ogra[i]=0.217578043021078


# In[140]:


x_major_locator = MultipleLocator(60)
y_major_locator = MultipleLocator(0.5)

fig, axes = plt.subplots(1,1)

axes.plot(tlist, evol0.expect[0],linewidth=1.5, label=r'$|\langle \psi(t)|\psi(0)\rangle|^2 $',color='silver')
axes.plot(tlist, evol0.expect[1],linewidth=1.5, label=r'$\langle \hat{O}\rangle$',color='blue')
axes.plot(tlist, Ocan,linewidth=1.5, linestyle='-.',label=r'$\rho_c(\beta_c)$',color='red')
axes.plot(tlist, Ogra,linewidth=1.5, label=r'$\rho(\beta,\alpha)$')
#plt.xlim(19,21)
#plt.ylim(0.3,0.8)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
axes.xaxis.set_major_locator(x_major_locator)
axes.yaxis.set_major_locator(y_major_locator)
#plt.xlabel('t',fontsize=20)
#plt.ylabel(r'$|\langle \psi(0) | \rho(t) |\psi(0) \rangle|^2$',fontsize=20)
axes.legend(loc=1)
#axes.legend(framealpha=0.5)
plt.rc('axes',linewidth=2)
plt.savefig('psi0uumdd.svg')


# In[ ]:




