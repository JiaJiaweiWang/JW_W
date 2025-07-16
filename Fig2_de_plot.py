
# In[1]:


import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from scipy import optimize as op
from matplotlib.ticker import MultipleLocator


# In[2]:


N=7               #system size
j=1            #spin size
d=int(2*j+1)      #


# In[3]:


idd=qeye(d)         #single site unit

ids=idd             #unite
for i in range(N-1):
    ids=tensor(ids,idd)

Em=np.zeros((d,d))

em=Qobj(Em)         #empty

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


# In[4]:


evolve1=qload('evolve1')
evolve2=qload('evolve2')
evolve3=qload('evolve3')
evolve4=qload('evolve4')
evolve5=qload('evolve5')
evolve6=qload('evolve6')


# In[5]:


Esample=qload('esample')
Vsample=qload('vsample')


# In[6]:


def st(i):
    state=basis(d,i)
    return(state)

rdd=idd
for i in range(N-3):
    rdd=tensor(rdd,idd)
py=(tensor(st(0),st(1))+tensor(st(1),st(2)))/(2**(1/2))
Py=tensor(py*py.dag(),rdd)
expy=[]
for i in range(9):
    expy.append(N*np.real((Vsample[i].dag()*Py*Vsample[i]).full()[0,0]))


# In[7]:


def f(x,A):
    return(np.exp(-A*x))


# In[8]:


tlist=np.linspace(0,40,400)


# In[157]:


drate=[]
emsf=[]
for i in range(9):
    As=op.curve_fit(f,tlist,evolve1[i].expect[0])[0]
    drate.append(As[0])
    emsf.append(f(tlist,As))


# In[158]:


evolveplot=evolve1


# In[266]:


x_major_locator = MultipleLocator(20)
y_major_locator = MultipleLocator(0.3)

fig, axes = plt.subplots(1,1)

axes.plot(tlist, evolveplot[0].expect[0],linewidth=2,color='silver')
axes.plot(tlist, evolveplot[1].expect[0],linewidth=2,color='silver')
axes.plot(tlist, evolveplot[2].expect[0],linewidth=2,color='silver',label= r'$scar states$')
axes.plot(tlist, evolveplot[3].expect[0],linewidth=1.5,color='deepskyblue')
axes.plot(tlist, evolveplot[4].expect[0],linewidth=1.5,color='deepskyblue')
axes.plot(tlist, evolveplot[5].expect[0],linewidth=1.5,color='deepskyblue')
axes.plot(tlist, evolveplot[6].expect[0],linewidth=1.5,color='deepskyblue')
axes.plot(tlist, evolveplot[7].expect[0],linewidth=1.5,color='deepskyblue')
axes.plot(tlist, evolveplot[8].expect[0],linewidth=1.5,color='deepskyblue',label= r'$thermal states$')
axes.plot(tlist, emsf[0],linewidth=2,linestyle='-.',color='dimgray')
axes.plot(tlist, emsf[1],linewidth=2,linestyle='-.',color='dimgray')
axes.plot(tlist, emsf[2],linewidth=2,linestyle='-.',color='dimgray',label= r'$fit e^{-\alpha_i t} for scar states$')
axes.plot(tlist, emsf[3],linewidth=1.6,linestyle='-.',color='crimson')
axes.plot(tlist, emsf[4],linewidth=1.6,linestyle='-.',color='crimson')
axes.plot(tlist, emsf[5],linewidth=1.6,linestyle='-.',color='crimson')
axes.plot(tlist, emsf[6],linewidth=1.6,linestyle='-.',color='crimson')
axes.plot(tlist, emsf[7],linewidth=1.6,linestyle='-.',color='crimson')
axes.plot(tlist, emsf[8],linewidth=1.6,linestyle='-.',color='crimson',label= r'$fit e^{-\alpha_i t} for thermal states$')
plt.xlim([0,40])
plt.ylim([0.4,1])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
axes.xaxis.set_major_locator(x_major_locator)
axes.yaxis.set_major_locator(y_major_locator)
plt.xlabel('t',fontsize=20)
plt.ylabel(r'$|\langle E_i | \rho(t) |E_i \rangle|^2$',fontsize=20)
axes.legend(loc=1)
axes.legend(framealpha=0.5,fontsize=15)
plt.rc('axes',linewidth=2)
plt.savefig('spin1EigenDecay3.svg')


# In[261]:


plt.plot(expy,drate,'o')


# In[268]:


def ff(x,A):
    return(A*x)
Ass=op.curve_fit(ff,expy,drate)[0]
xn=[0,1.5]
dratefit=ff(xn,Ass)

x_major_locator = MultipleLocator(1)
y_major_locator = MultipleLocator(0.01)

fig, axes = plt.subplots(1,1)
plt.plot(xn,dratefit,linewidth=4,linestyle=':',label= r'$\alpha=\kappa_{fit}\cdot\langle \hat{N}\rangle$',color='blue')
plt.plot(expy,drate,'o',markersize=8,label= r'$(\langle \hat{N}\rangle_i,\alpha_{i})$',color='black')
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlim([0,1.5])
plt.ylim([0,0.021])
axes.xaxis.set_major_locator(x_major_locator)
axes.yaxis.set_major_locator(y_major_locator)
axes.ticklabel_format(style='sci',scilimits=(-1,2),axis='y')
axes.get_yaxis().get_offset_text().set(va='bottom',ha='left')
axes.yaxis.get_offset_text().set_fontsize(25)
plt.xlabel(r'$\langle \hat{N} \rangle$',fontsize=25)
plt.ylabel(r'$\alpha$',fontsize=25)
axes.legend(framealpha=0.5,fontsize=20,loc=2)
plt.rc('axes',linewidth=3.5)
plt.savefig('alphafit.svg')


# In[245]:
#ass=[]
ass.append(Ass[0])  # this needs to be repeated for different coefficient c, then 'ass' can be utilized in the calculation below.



cinverse=np.zeros((6,1))
for i in range(6):
    cinverse[i,0]=1/((i+1)*200)

cin=[0.0007,0.0052]
kes=np.zeros((2,1))
for i in range(2):
    kes[i,0]=2*((2*j)**(1/2))*cin[i]


# In[251]:


x_major_locator = MultipleLocator(0.002)
y_major_locator = MultipleLocator(0.005)

fig, axes = plt.subplots(1,1)
plt.plot(cin,kes,linewidth=3,label=r'$2\sqrt{2j}\cdot 1/c$',color='green')
plt.plot(cinverse,ass,'o',markersize=8,label=r'$\kappa_{fit}$',color='blue')
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlim([cin[0],cin[1]])
plt.ylim([kes[0,0],kes[1,0]])
axes.xaxis.set_major_locator(x_major_locator)
axes.yaxis.set_major_locator(y_major_locator)
axes.ticklabel_format(style='sci',scilimits=(-1,2),axis='y')
axes.get_yaxis().get_offset_text().set(va='bottom',ha='left')
axes.yaxis.get_offset_text().set_fontsize(25)
axes.ticklabel_format(style='sci',scilimits=(-1,2),axis='x')
axes.get_xaxis().get_offset_text().set(va='bottom',ha='left')
axes.xaxis.get_offset_text().set_fontsize(25)
plt.xlabel(r'$1/c$',fontsize=25)
plt.ylabel(r'$\kappa$',fontsize=25)
axes.legend(framealpha=1,loc=2,fontsize=20)
plt.rc('axes',linewidth=3)
plt.savefig('kappacompare.svg')



