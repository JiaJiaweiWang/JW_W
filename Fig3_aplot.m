load('Oexpectation','Oexp')
load('en9spin1','En')
load('dim9spin1','num')
Ene=-7.5:0.5:7.5;
E=zeros(num,1);
for i=1:num
    E(i,1)=real(En(i,i));
end
Ocano=zeros(size(Ene,2),1);
for i=1:size(Ene,2)
    betac=BetaC(Ene(i));
    ocu=0;
    ocd=0;
    for k=1:num
        ocu=ocu+exp(-betac*real(En(k,k)))*Oexp(k,1);
        ocd=ocd+exp(-betac*real(En(k,k)));
    end
    Ocano(i,1)=ocu/ocd;
end
fit1 = polyfit(Ene, Ocano, 10); % p 返回十次多项式的系数
Enne=min(E):0.1:max(E);
Ocanofit=polyval(fit1,Enne);
Ocdeviat=zeros(num,1);
for i=1:num
    Ocdeviat(i,1)=real(Oexp(i,1)-polyval(fit1,E(i,1)));
end
scatter(E,Oexp,40,Ocdeviat,'filled','MarkerEdgeColor', 'k')
c = colorbar('Location', 'westoutside'); % 显示颜色条
caxis([min(Ocdeviat) max(Ocdeviat)]);
c.Ticks = linspace(-0.2 , 0.3,3); % 设置颜色条刻度
colormap(turbo)

hold on
plot(Enne,Ocanofit,'LineWidth',1.2,'color','b','LineStyle','-.')
xlim([-7.7,7.7])
ylim([-0.1,0.6])
xticks(-5:5:5)
yticks(-0.1:0.3:0.5)
box on
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 300, 500, 400]);


function [beta]=BetaC(E)
load('en9spin1','En');
num=size(En,1);
beta=fzero(@EcEnn,[-10,10]);
function [dE]=EcEnn(beta1)
    fu=0;
    fd=0;
    for l=1:num
        fu=fu+exp(-beta1*En(l,l))*En(l,l);
        fd=fd+exp(-beta1*En(l,l));
    end
    Ec=fu/fd;
    dE=Ec-E;
end
end