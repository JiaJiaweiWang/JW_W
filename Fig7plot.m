global num En Np Eig pye
% %基本算符定义
j=1;
d=round(4*j+1);
N=6;
ids=speye(d,d);
rdd1=speye(d^(N-1),d^(N-1));
rdd2=speye(d^(N-2),d^(N-2));
idd=speye(d^N,d^N);

% %Py的定义
R=zeros(d,d);
L=zeros(d,d);
for i=1:2*j
    R(i,i)=1;
    L(d-i+1,d-i+1)=1;
end
r=sparse(R);
l=sparse(L);
h0=zeros(d,d);
for i=1:2*j
    m=j-i;
    hm=(1/2)*(((j-m)*(j+m+1))^(1/2));
    h0(i,i+1)=hm;
    h0(i+1,i)=hm;
    h0(d-i+1,d-i)=hm;
    h0(d-i,d-i+1)=hm;
end
h=sparse(h0);
r1=chainOp(r,1,N);
l2=chainOp(l,2,N);
h1=chainOp(h,1,N);
h2=chainOp(h,2,N);
Py=(h1+h2)*r1*l2*(h1+h2);

Np=expectationO(Py);
% save('Npy','Np')
% scatter(E,Np)

load('dim12spin1','num');
load('en12spin1','En');
% load('Npy','Np');
alp=zeros(num,1);
bet=zeros(num,1);
% % % % % %利用Py_i与E_i求出相应的beta，alpha
x0=[0,0];
for k=1:num
    Eig=En(k,k);
    pye=Np(k,1);
    x=fsolve(@disE_Pyg,x0);
    bet(k,1)=x(1);
    alp(k,1)=x(2);
end

betc=zeros(num,1);
x1=0;
for k=1:num
    Eig=En(k,k);
    betc(k,1)=fsolve(@disEonly,x1);
end

save('beta','bet')
save('alpha','alp')
save('betaC','betc')

% load('Npy.mat','Np')
% load('en12spin1.mat','En')
% load('dim12spin1.mat','num')
% load('beta.mat','bet')
% load('alpha.mat','alp')
% load('betaC.mat','betc')
% load('Snum','snum')
% load('Nsnum','nsnum')


% rdd1=speye(5^5,5^5);
% rdd2=speye(5^4,5^4);
% sz=sparse([1,0,0,0,0;0,0,0,0,0;0,0,-1,0,0;0,0,0,1,0;0,0,0,0,1]);
% Ob=kron(sz,rdd1);
% expr=expectationO(Ob);


expr=Np;
expg=zeros(num,1);
expc=zeros(num,1);
E=zeros(num,1);
for i=1:num
    expg(i,1)=expG(bet(i,1),alp(i,1),expr);
    expc(i,1)=expC(betc(i,1),expr);
    E(i,1)=En(i,i);
end

Eexpc=[E,expc];
Eexpcs=sortrows(Eexpc,'ComparisonMethod','real');

% scatter(E(snum),expr(snum))
% hold on
% scatter(E(nsnum),expr(nsnum))

po = [];
method = [];
radius = [];
N = [];
n = [];
ms = 25;
scatplot(real(E(nsnum)),real(expr(nsnum)),method,radius,N,n,po,ms);
hold on 
scatter(real(E(snum)),real(expr(snum)),40,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[0 0 1],'LineWidth',2) 
hold on 
plot(Eexpcs(:,1),Eexpcs(:,2),'LineWidth',2,'color',[0 0 0.8],'LineStyle','-.')
ylim([0.05,0.25])
xlim([-6.5,6.5])
xticks(-6:4:6)
yticks(0.1:0.1:0.2)
% % xlabel('E','FontName','Arial')
% % ylabel('log_{10}|\langle E_i|\psi(0)\rangle|^2','FontName','Arial')
box on
colorbar off
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 500, 450, 400]);

expgr=zeros(num,1);
expcr=zeros(num,1);
for i=1:num
    expgr(i,1)=expg(i,1)-expr(i,1);
    expcr(i,1)=expc(i,1)-expr(i,1);
end

scatter(E,expcr)
hold on
scatter(E,expgr,'filled')

scatter(E(nsnum),expcr(nsnum),'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',1.5)
hold on
scatter(E(snum),expcr(snum),'MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',1.5)
hold on
scatter(E(nsnum),expgr(nsnum),'MarkerEdgeColor',[0.05,0.21,0.92],'MarkerFaceColor',[0.31,0.78,0.96],'LineWidth',1)
hold on
scatter(E(snum),expgr(snum),'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[0.31,0.78,0.96],'LineWidth',1.5)
ylim([-0.16,0.12])
xlim([-6.5,6.5])
xticks(-6:4:6)
yticks(-0.1:0.1:0.1)
% xlabel('E','FontName','Arial')
% ylabel('log_{10}|\langle E_i|\psi(0)\rangle|^2','FontName','Arial')
box on
% colorbar off
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 500, 450, 400]);

function [expo]=expectationO(O)
load('vn12spin1.mat','Vn')
load('dim12spin1.mat','num')
expo=zeros(num,1);
for i=1:num
    expo(i,1)=Vn(:,i)'*O*Vn(:,i);
end
end

function [expO]=expG(beta,alpha,expr)
load('dim12spin1.mat','num')
load('en12spin1.mat','En')
load('Npy.mat','Np')
up=0;
down=0;
for i=1:num
    expon=beta*En(i,i)-alpha*Np(i,1);
    up=up+exp(-expon)*expr(i,1);
    down=down+exp(-expon);
end
expO=up/down;
end

function [expO]=expC(beta,expr)
load('dim12spin1.mat','num')
load('en12spin1.mat','En')
up=0;
down=0;
for i=1:num
    expon=beta*En(i,i);
    up=up+exp(-expon)*expr(i,1);
    down=down+exp(-expon);
end
expO=up/down;
end

% % % % 求解beta，alpha用到的，给出密度矩阵对能量和Py之预言的偏离的方程
function dEdPy=disE_Pyg(bet_alp)
global num En Np Eig pye
Efu=0;
Nfu=0;
fd=0;
for i=1:num
    expon=bet_alp(1)*En(i,i)-bet_alp(2)*Np(i,1);
    Efu=Efu+exp(-expon)*En(i,i);
    Nfu=Nfu+exp(-expon)*Np(i,1);
    fd=fd+exp(-expon);
end
Eg=Efu/fd;
Pyg=Nfu/fd;
dEdPy(1)=real(Eg-Eig);
dEdPy(2)=real(Pyg-pye);
end

% % % % 求解betac用到的，给出密度矩阵对能量和Py之预言的偏离的方程
function dE=disEonly(betonly)
global num En Eig
Efu=0;
fd=0;
for i=1:num
    expon=betonly*En(i,i);
    Efu=Efu+exp(-expon)*En(i,i);
    fd=fd+exp(-expon);
end
Ec=Efu/fd;
dE=real(Ec-Eig);
end