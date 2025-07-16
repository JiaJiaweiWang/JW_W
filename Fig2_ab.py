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


# In[7]:


Hp=Pb*Pt*H*Pt*Pb

Ei=Hp.eigenstates()[0]
Vi=Hp.eigenstates()[1]

En=[]
Vn=[]
for i in range(d**N):
    if abs(Ei[i])>10**(-12):
        En.append(Ei[i])
        Vn.append(Vi[i])


# In[8]:


Py=(tensor(st(0),st(1))+tensor(st(1),st(2)))/(2**(1/2))
Ny=Py*Py.dag()
for i in range(N-2):
    Ny=tensor(Ny,idd)


# In[9]:


psi0=tensor(st(0),st(0),st(0),st(0),st(0),st(0),st(0))
ovl=np.zeros(len(En))
Ne=np.zeros(len(En))
for i in range(len(En)):
    if abs((psi0.dag()*Vn[i]).full()[0,0])!=0:
        ovl[i]=np.log(abs((psi0.dag()*Vn[i]).full()[0,0]))
        Ne[i]=7*((Vn[i].dag()*Ny*Vn[i]).full()[0,0]).real
    else:
        ovl[i]=-10
        Ne[i]=7*((Vn[i].dag()*Ny*Vn[i]).full()[0,0]).real


# In[10]:


Es=[]
Vs=[]
Ens=[]
Vns=[]
ovls=[]
ovlns=[]
Ns=[]
Nns=[]
for i in range(len(En)):
    e=En[i]
    ov=-0.1*e**2-1.6
    if ovl[i]>ov:
        Es.append(En[i])
        Vs.append(Vn[i])
        ovls.append(ovl[i])
        Ns.append(Ne[i])
    else:
        Ens.append(En[i])
        Vns.append(Vn[i])
        ovlns.append(ovl[i])
        Nns.append(Ne[i])


# In[12]:


Es


# In[19]:


es=Es[6]
vs=Vs[6]
ovs=ovls[6]
ns=Ns[6]

ens=[]
vns=[]
ovns=[]
nns=[]

for i in range(len(Ens)):
    if(Ens[i]<-0.7 and Ens[i]>-1.2 and ovlns[i]>-7):
        ens.append(Ens[i])
        vns.append(Vns[i])
        ovns.append(ovlns[i])
        nns.append(Nns[i])


# In[20]:


ens


# In[25]:


emax=[-0.7,-0.7]
emin=[-1.2,-1.2]
oedge=[-7,0]

x_major_locator = MultipleLocator(4)
y_major_locator = MultipleLocator(4)

fig, axes = plt.subplots(1,1)
axes.scatter(Es,ovls, facecolors='blue', edgecolors='darkgray', s=100)
axes.scatter(Ens,ovlns, facecolors='darkgray', edgecolors='gray', s=100)
axes.plot(emax,oedge,linestyle=':',color='red',linewidth=3)
axes.plot(emin,oedge,linestyle=':',color='red',linewidth=3)
axes.scatter(ens,ovns, facecolors='none', edgecolors='black', s=200, linewidths=3)
axes.scatter(es,ovs, facecolors='none', edgecolors='red', s=200, linewidths=3)
plt.ylim([-7,0])
#axes.xaxis.tick_top()
#axes.xaxis.set_label_position('top')
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
axes.xaxis.set_major_locator(x_major_locator)
axes.yaxis.set_major_locator(y_major_locator)
plt.xlabel('E',fontsize=20)
plt.ylabel('overlap',fontsize=20)
plt.rc('axes',linewidth=2)
plt.savefig('ovlp.svg')


# In[23]:


vt=[]
vt.append(vs)
for i in range(8):
    vt.append(vns[i])
    
qsave(vt,'vt')


# In[24]:


emax=[-0.7,-0.7]
emin=[-1.2,-1.2]
oedge=[0,2]

x_major_locator = MultipleLocator(4)
y_major_locator = MultipleLocator(1)

fig, axes = plt.subplots(1,1)
axes.scatter(Es, Ns, facecolors='blue', edgecolors='darkblue', s=100)
axes.scatter(Ens, Nns, facecolors='darkgray', edgecolors='gray', s=100)
axes.plot(emax,oedge,linestyle=':',color='red',linewidth=3)
axes.plot(emin,oedge,linestyle=':',color='red',linewidth=3)
axes.scatter(ens, nns, facecolors='none', edgecolors='black', s=200, linewidths=3)
axes.scatter(es, ns, facecolors='none', edgecolors='red', s=200, linewidths=3)
plt.ylim([0,2])
#axes.xaxis.tick_top()
#axes.xaxis.set_label_position('top')
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
axes.xaxis.set_major_locator(x_major_locator)
axes.yaxis.set_major_locator(y_major_locator)
plt.xlabel('E',fontsize=20)
plt.ylabel('overlap',fontsize=20)
plt.rc('axes',linewidth=2)
plt.savefig('Nexp.svg')





