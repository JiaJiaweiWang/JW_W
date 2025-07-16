
# In[1]:


import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from scipy import optimize as op
import matplotlib.pyplot as plt


# In[2]:


N=7               #链长
j=1            #自旋大小
d=int(2*j+1)      #自旋矩阵维度


# In[34]:


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


# In[35]:


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


# In[45]:


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


# In[46]:


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


# In[38]:


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


# In[40]:


Hp=Pb*Pt*H*Pt*Pb


# In[39]:


Ei=Hp.eigenstates()[0]
Vi=Hp.eigenstates()[1]

En=[]
Vn=[]
for i in range(d**N):
    if abs(Ei[i])>10**(-12):
        En.append(Ei[i])
        Vn.append(Vi[i])


# In[11]:


psi0=tensor(st(0),st(0),st(0),st(0),st(0),st(0),st(0))
ovl=np.zeros(len(En))
for i in range(len(En)):
    if abs((psi0.dag()*Vn[i]).full()[0,0])!=0:
        ovl[i]=np.log(abs((psi0.dag()*Vn[i]).full()[0,0]))
    else:
        ovl[i]=-10


# In[12]:


Es=[]
Vs=[]
Ens=[]
Vns=[]
ovls=[]
ovlns=[]
for i in range(len(En)):
    e=En[i]
    ov=-0.1*e**2-1.6
    if ovl[i]>ov:
        Es.append(En[i])
        Vs.append(Vn[i])
        ovls.append(ovl[i])
    else:
        Ens.append(En[i])
        Vns.append(Vn[i])
        ovlns.append(ovl[i])


# In[17]:


Vsample=[]
Esample=[]
for i in range(3):
    Vsample.append(Vs[2*i+1])
    Esample.append(Es[2*i+1])
for i in range(6):
    Vsample.append(Vns[8*i])
    Esample.append(Ens[8*i])


# In[20]:


qsave(Vsample,'vsample')
qsave(Esample,'esample')


# In[6]:


Vsample=qload('vsample')
Esample=qload('esample')


# In[28]:


tlist=np.linspace(0,40,400)
trajn=500


# In[47]:


evolve1=[]
cn=200
Lindc=Lind(cn)
psi0=Vsample[0]
rho0=psi0*psi0.dag()
emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn)
opts = Options(rhs_reuse=True)
evolve1.append(emc)

for i in range(8):
    psi0=Vsample[i+1]
    rho0=psi0*psi0.dag()
    emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn, options=opts)
    evolve1.append(emc)


# In[48]:


evolve2=[]
cn=400
Lindc=Lind(cn)
psi0=Vsample[0]
rho0=psi0*psi0.dag()
emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn)
opts = Options(rhs_reuse=True)
evolve2.append(emc)

for i in range(8):
    psi0=Vsample[i+1]
    rho0=psi0*psi0.dag()
    emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn, options=opts)
    evolve2.append(emc)


evolve3=[]
cn=600
Lindc=Lind(cn)
psi0=Vsample[0]
rho0=psi0*psi0.dag()
emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn)
opts = Options(rhs_reuse=True)
evolve3.append(emc)

for i in range(8):
    psi0=Vsample[i+1]
    rho0=psi0*psi0.dag()
    emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn, options=opts)
    evolve3.append(emc)



evolve4=[]
cn=800
Lindc=Lind(cn)
psi0=Vsample[0]
rho0=psi0*psi0.dag()
emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn)
opts = Options(rhs_reuse=True)
evolve4.append(emc)

for i in range(8):
    psi0=Vsample[i+1]
    rho0=psi0*psi0.dag()
    emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn, options=opts)
    evolve4.append(emc)



evolve5=[]
cn=1000
Lindc=Lind(cn)
psi0=Vsample[0]
rho0=psi0*psi0.dag()
emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn)
opts = Options(rhs_reuse=True)
evolve5.append(emc)

for i in range(8):
    psi0=Vsample[i+1]
    rho0=psi0*psi0.dag()
    emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn, options=opts)
    evolve5.append(emc)


evolve6=[]
cn=1200
Lindc=Lind(cn)
psi0=Vsample[0]
rho0=psi0*psi0.dag()
emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn)
opts = Options(rhs_reuse=True)
evolve6.append(emc)

for i in range(8):
    psi0=Vsample[i+1]
    rho0=psi0*psi0.dag()
    emc=mcsolve(H0, psi0, tlist,Lindc, [rho0], ntraj=trajn, options=opts)
    evolve6.append(emc)


qsave(evolve1,'evolve1')
qsave(evolve2,'evolve2')
qsave(evolve3,'evolve3')
qsave(evolve4,'evolve4')
qsave(evolve5,'evolve5')
qsave(evolve6,'evolve6')

