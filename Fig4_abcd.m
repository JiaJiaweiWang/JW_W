global num En py_tunning Eig pye
load('dim9spin1','num');
load('en9spin1','En');
load('Py_tunning','py_tunning');
load('ens','Ens')
load('es','Es')
load('vs','Vs')
load('vns','Vns')
load('pyns','expns')
load('pys','exps')
load('Vsubs','vsubs')
load('Vsubns','vsubns')
load('vsub9spin1','Vsub')
snum=size(Es,1);
nsnum=size(Ens,1);
alps=zeros(snum,1);
bets=zeros(snum,1);
alpns=zeros(nsnum,1);
betns=zeros(nsnum,1);
% % % % %利用Py_i与E_i求出相应的beta，alpha
x0=[0,0];
for k=1:snum
    Eig=Es(k,1);
    pye=exps(k,1);
    x=fsolve(@disE_Pyg,x0);
    bets(k,1)=x(1);
    alps(k,1)=x(2);
end
for k=1:nsnum
    Eig=Ens(k,1);
    pye=expns(k,1);
    x=fsolve(@disE_Pyg,x0);
    betns(k,1)=x(1);
    alpns(k,1)=x(2);
end


% % % % 利用前述的alpha和beta求出一些可观测量的期望值，并与本征态的结果相比较
d=9;
id3=speye(3,3);
rdd1=speye(3^(d-1),3^(d-1));
rdd2=speye(3^(d-2),3^(d-2));
idd=speye(3^d,3^d);
s1=sparse([1;0;0]);
s2=sparse([0;1;0]);
s3=sparse([0;0;1]);
sz=sparse([1,0,0;0,0,0;0,0,-1]);
ob=s1*s1'+s2*s2';
Ob=kron(ob,id3);
Obw=kron(ob,rdd1);
expsg=zeros(snum,1);
expnsg=zeros(nsnum,1);
expsc=zeros(snum,1);
expnsc=zeros(nsnum,1);
expsr=zeros(snum,1);
expnsr=zeros(nsnum,1);
for i=1:snum
    if(i==1)
        Ee=Es(i,1)-0.01;
    elseif(i==2)
        Ee=Es(i,1)+0.01;
    else
        Ee=Es(i,1);
    end
    expsr(i,1)=real(Vs(:,i)'*Obw*Vs(:,i));
    expsg(i,1)=trace(Ob*Vsub_grand(bets(i,1),alps(i,1)));
    expsc(i,1)=trace(Ob*Vsub_grand(BetaCPy(Ee,0),0));
end
for i=1:nsnum
    expnsr(i,1)=real(Vns(:,i)'*Obw*Vns(:,i));
    expnsg(i,1)=trace(Ob*Vsub_grand(betns(i,1),alpns(i,1)));
    expnsc(i,1)=trace(Ob*Vsub_grand(BetaCPy(Ens(i,1),0),0));
end
dsgr=zeros(snum,1);
dscr=zeros(snum,1);
dnsgr=zeros(nsnum,1);
dnscr=zeros(nsnum,1);
for i=1:snum
    dsgr(i,1)=expsg(i,1)-expsr(i,1);
    dscr(i,1)=expsc(i,1)-expsr(i,1);
end
for i=1:nsnum
    dnsgr(i,1)=expnsg(i,1)-expnsr(i,1);
    dnscr(i,1)=expnsc(i,1)-expnsr(i,1);
end

scatter(Ens,dnscr,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',1.5)
hold on
scatter(Es,dscr,'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1.5)
hold on
scatter(Ens,dnsgr,'MarkerEdgeColor',[0.05,0.21,0.92],'MarkerFaceColor',[0.31,0.78,0.96],'LineWidth',1)
hold on
scatter(Es,dsgr,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[0.31,0.78,0.96],'LineWidth',1.5)
ylim([-0.04,0.08])
xlim([-8,8])
xticks(-8:4:8)
yticks(-0.04:0.04:0.08)
% xlabel('E','FontName','Arial')
% ylabel('log_{10}|\langle E_i|\psi(0)\rangle|^2','FontName','Arial')
box on
% colorbar off
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 500, 450, 400]);

vsubs=zeros(9,9,snum);
vsubns=zeros(9,9,nsnum);
ks=0;
kns=0;
for i=1:num
    of=-0.004*E(i,1)^2+0.27;
    if of<ooo(i,1)
        ks=ks+1;
        vsubs(:,:,ks)=Vsub(:,:,i);
    else
        kns=kns+1;
        vsubns(:,:,kns)=Vsub(:,:,i);
    end
end


% % % % 预期的与真实的本征态的局域约化密度矩阵的距离
p=1;
dsg=zeros(snum,1);
dsc=zeros(snum,1);
dnsg=zeros(nsnum,1);
dnsc=zeros(nsnum,1);
for i=1:snum
    if(i==1)
        Ee=Es(i,1)-0.01;
    elseif(i==2)
        Ee=Es(i,1)+0.01;
    else
        Ee=Es(i,1);
    end    
    vsubg=Vsub_grand(bets(i,1),alps(i,1));
    vsubc=Vsub_grand(BetaCPy(Ee,0),0);
    dsg(i,1)=sdistance(vsubg,vsubs(:,:,i),p);
    dsc(i,1)=sdistance(vsubc,vsubs(:,:,i),p);
end
for i=1:nsnum
    vsubg=Vsub_grand(betns(i,1),alpns(i,1));
    vsubc=Vsub_grand(BetaCPy(Ens(i,1),0),0);
    dnsg(i,1)=sdistance(vsubg,vsubns(:,:,i),p);
    dnsc(i,1)=sdistance(vsubc,vsubns(:,:,i),p);
end
scatter(Ens,dnsc,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',1.5)
hold on
scatter(Es,dsc,'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1.5)
hold on
scatter(Ens,dnsg,'MarkerEdgeColor',[0.05,0.21,0.92],'MarkerFaceColor',[0.31,0.78,0.96],'LineWidth',1)
hold on
scatter(Es,dsg,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[0.31,0.78,0.96],'LineWidth',1.5)
ylim([0,1.2])
xlim([-8,8])
xticks(-8:4:8)
yticks(0:0.4:1.2)
% xlabel('E','FontName','Arial')
% ylabel('log_{10}|\langle E_i|\psi(0)\rangle|^2','FontName','Arial')
box on
% colorbar off
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 500, 450, 400]);


% % solving for beta and alpha with energy E and quasi-particle number N
function dEdPy=disE_Pyg(bet_alp)
global num En py_tunning Eig pye
Efu=0;
Nfu=0;
fd=0;
for i=1:num
    expon=bet_alp(1)*En(i,i)-bet_alp(2)*py_tunning(i,1);
    Efu=Efu+exp(-expon)*En(i,i);
    Nfu=Nfu+exp(-expon)*py_tunning(i,1);
    fd=fd+exp(-expon);
end
Eg=Efu/fd;
Pyg=Nfu/fd;
dEdPy(1)=real(Eg-Eig);
dEdPy(2)=real(Pyg-pye);
end

% grand canonical ensemble of subsystem under thermodynamical parameter beta and alpha
function [vsubg]=Vsub_grand(beta,alpha)
load('vsub9spin1','Vsub');
load('dim9spin1','num');
load('en9spin1','En');
load('Py_tunning','py_tunning');
vsubgu=zeros(9,9);
vsubgd=0;
for i=1:num
    expon=beta*En(i,i)-alpha*py_tunning(i,1);
    vsubgu=vsubgu+exp(-expon)*Vsub(:,:,i);
    vsubgd=vsubgd+exp(-expon);
end
vsubg=vsubgu/vsubgd;
end

function [beta]=BetaCPy(E,alpha)
load('en9spin1','En');
load('N','N')
num=size(En,1);
beta=fzero(@EcEnn,[-20,20]);
function [dE]=EcEnn(beta1)
    fu=0;
    fd=0;
    for l=1:num
        fu=fu+exp(-beta1*En(l,l)+alpha*N(l,1))*En(l,l);
        fd=fd+exp(-beta1*En(l,l)+alpha*N(l,1));
    end
    Ec=fu/fd;
    dE=Ec-E;
end
end