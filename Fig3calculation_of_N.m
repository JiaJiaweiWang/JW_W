load('dim9spin1','num')
load('vn9spin1','Vn');
load('en9spin1','En');
load('vsubhalf9spin1','Vsub');
d=3;
N=9;
id4=speye(d,d);
rdd1=speye(d^(N-1),d^(N-1));
rdd2=speye(d^(N-2),d^(N-2));
idd=speye(d^N,d^N);
st1=sparse([1;0;0]);
st2=sparse([0;1;0]);
st3=sparse([0;0;1]);
inia=kron(st1,kron(st1,kron(st1,kron(st1,kron(st1,kron(st1,kron(st1,kron(st1,st1))))))));
ove=zeros(num,1);
E=zeros(num,1);
for i=1:num
    ove(i,1)=log10(abs(inia'*Vn(:,i))^2);
    E(i,1)=real(En(i,i));
end
% scatter(E,ove,'filled');
% xlim([-9,9])
% ylim([-11,0])
% xf=-9:0.01:9;
% yf=-0.0664*xf.^2-1.25;
% scatter(E,ove)
% hold on
% plot(xf,yf)
snum=18;
nsnum=num-snum;
Es=zeros(snum,1);
Ens=zeros(nsnum,1);
Vs=zeros(3^9,snum);
Vns=zeros(3^9,nsnum);
Os=zeros(snum,1);
Ons=zeros(nsnum,1);
vsubs=zeros(3^4,3^4,snum);
vsubns=zeros(3^4,3^4,nsnum);
ks=0;
kns=0;
for i=1:num
    of=-0.0664*E(i,1)^2-1.25;
    if of<ove(i,1)
        ks=ks+1;
        Es(ks,1)=E(i,1);
        Os(ks,1)=ove(i,1);
        Vs(:,ks)=Vn(:,i);
        vsubs(:,:,ks)=Vsub(:,:,i);
    else
        kns=kns+1;
        Ens(kns,1)=E(i,1);
        Ons(kns,1)=ove(i,1);
        Vns(:,kns)=Vn(:,i);
        vsubns(:,:,kns)=Vsub(:,:,i);
    end
end
% scatter(Es,Os)
% hold on
% scatter(Ens,Ons)
save('vs','Vs');
save('vns','Vns');
save('es','Es');
save('ens','Ens');
save('Vsubhalfs','vsubs');
save('Vsubhalfns','vsubns');
 
ob=(kron(st1,st2)+kron(st2,st3))/(2^(1/2));
% ob=kron(st1*st1',st2*st2')+kron(st2*st2',st3*st3');
Ob=kron(ob*ob',rdd2);
exps=zeros(snum,1);
expns=zeros(nsnum,1);
for i=1:snum
    exps(i,1)=Vs(:,i)'*Ob*Vs(:,i);
end
for i=1:nsnum
    expns(i,1)=real(Vns(:,i)'*Ob*Vns(:,i));
end
% scatter(Ens,expns)
% hold on
% scatter(Es,exps)
save('pys','exps');
save('pyns','expns');

Ene=-7.5:0.5:7.5;
pymean=zeros(size(Ene,2),1);
for i=1:size(Ene,2)
    pymean(i,1)=PrMean(Ene(i),Ob);
end
% scatter(Ens,expns) 
% hold on 
% scatter(Es,exps) 
% hold on 
% plot(Ene,pymean)
% 
% % % % % % 孤立系统中改进ETH的尝试
py_tunning=zeros(num,1);
psiy=(kron(st1,st2)+kron(st2,st3))/(2^(1/2));
Py=kron(psiy*psiy',rdd2);
for i=1:num
    py_tunning(i,1)=real(Vn(:,i)'*Py*Vn(:,i));
end
save('Py_tunning','py_tunning')
% 
% pymean=zeros(size(Ene,2),5);
% for k=1:5
%     for i=1:size(Ene,2)
%         pymean(i,k)=PrMeanPy(Ene(i),Ob,-100+k*40);
%     end
% end
% scatter(Ens,expns)
% hold on
% scatter(Es,exps)
% hold on
% for k=1:5
%     plot(Ene,pymean(:,k))
%     hold on
% end
% 
% save('vscar','Vs')
% save('vnscar','Vns')
% save('escar','Es')
% save('enscar','Ens')
% save('pyscarexpect','exps')
% save('pynscarexpect','expns')
% 
% % % % Py在Hscar的本征态下的展开
% Eleft=zeros(num,num);
% Eright=zeros(num,num);
% Pym=zeros(num,num);
% for i=1:num
%     for k=1:num
%         Eleft(i,k)=En(i,i);
%         Eright(k,i)=En(i,i);
%         Pym(i,k)=Vn(:,i)'*Py*Vn(:,k);
%     end
% end
% scatter3(real(Eleft(:)),real(Eright(:)),imag(Pym(:)))
% 
% 
% pymeanR=zeros(size(Ene,2),1);
% for i=1:size(Ene,2)
%     pymeanR(i,1)=exps(1,1)+3*(pymean(i,1)-exps(1,1));
% end
% scatter(Ens,expns)
% hold on
% scatter(Es,exps)
% hold on
% plot(Ene,pymeanR)
% 
% 
% rati=zeros(snum-2,1);
% for i=1:snum-2
%     rati(i,1)=(pymean(i,1)-exps(1,1))/(exps(i+2,1)-exps(1,1));
% end
% 
% scatter(Esm,rati)