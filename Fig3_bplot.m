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
% % % % %求解
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


save('betas','bets')
save('betans','betns')
save('alphas','alps')
save('alphans','alpns')


% scatter(alps,bets)
% hold on
% scatter(alpns,betns)

id3=speye(3,3);
s1=sparse([1;0;0]);
s2=sparse([0;1;0]);
s3=sparse([0;0;1]);
sz=sparse([1,0,0;0,0,0;0,0,-1]); 
Ob=kron(sz,sz);
Oexp=zeros(num,1);
for i=1:num
    Oexp(i,1)=trace(Vsub(:,:,i)*Ob);
end

Oexps=zeros(18,1);
Oexpns=zeros(630,1);
for i=1:18
    Oexps(i,1)=trace(vsubs(:,:,i)*Ob);
end
for i=1:630
    Oexpns(i,1)=trace(vsubns(:,:,i)*Ob);
end
save('Osexpectation','Oexps')
save('Onsexpectation','Oexpns')
save('Oexpectation','Oexp')



% % % % % 计算函数值
Z=zeros(105,141);
for i=1:105
    for j=1:141
        up=0;
        down=0;
        for k=1:num
            expon=Bet(i,j)*(real(En(k,k))-Alph(i,j)*real(py_tunning(k,1)));
            up=up+exp(-expon)*real(Oexp(k,1));
            down=down+exp(-expon);
        end
        Z(i,j)=up/down;
    end
end

E=zeros(num,1);
for i=1:num
    E(i,1)=real(En(i,i));
end
scatter3(E,py_tunning,Oexp)


ft=fit([E,py_tunning],Oexp,'poly23');

Odeviat=real(Oexp-ft(E,py_tunning));

% 绘制结果
% plot(ft, [E, py_tunning], Oexp);
zzeros=zeros(18,1);
zzerons=zeros(630,1);
for i=1:18
    zzeros(i,1)=-0.2;
end
for i=1:630
    zzerons(i,1)=-0.2;
end
[Eg,Ng]=meshgrid(linspace(min(E), max(E), 20), linspace(min(py_tunning), max(py_tunning), 20));
Og=ft(Eg,Ng);
% 
% hold on
scatter3(Es,exps,zzeros,15,'k','filled', 'MarkerEdgeColor', 'r','MarkerFaceAlpha', 0.2,'MarkerEdgeAlpha', 0.8)
hold on
scatter3(Ens,expns,zzerons,15,'k','filled', 'MarkerEdgeColor', 'k','MarkerFaceAlpha', 0.2,'MarkerEdgeAlpha', 0.2)
hold on
scatter3(E,py_tunning,Oexp,40,Odeviat,'filled', 'MarkerEdgeColor', 'k')
xlim([min(E), max(E)])
ylim([min(py_tunning), max(py_tunning)])
c = colorbar('Location', 'eastoutside'); % 显示颜色条
caxis([min(Odeviat) max(Odeviat)]);
c.Ticks = linspace(-0.07 , 0.07,3); % 设置颜色条刻度
colormap(turbo)
% colorbar off
hold on
surf(Eg, Ng, Og, 'FaceAlpha', 0.1, 'EdgeColor', 'b','FaceColor','black','EdgeAlpha', 0.5)
%  'MeshStyle', 'Row'
xticks(-5:5:5)
yticks(0.1:0.1:0.2)
zticks(-0.1:0.3:0.5)
box on
set(gca,'linewidth',1.1)
set(gca,'FontName','Times New Roman','FontSize',18)
set(gcf, 'Position', [500, 300, 500, 400]);
% % % % % 绘制三维曲面图
% mesh(Alph,Bet,Z);
% hold on
% scatter3(alpns,betns,Oexpns,30,'filled')
% hold on
% scatter3(alps,bets,Oexps,30,'filled')




% 添加标签和标题
% xlabel('X轴');
% ylabel('Y轴');
% zlabel('Z轴');
% title('函数 Z = sin(X) * cos(Y) 的三维曲面图');


% % % % % 根据E，N求得alpha与beta
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