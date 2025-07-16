% j=1;
% d=round(4*j+1);
% No=12;
% N=round(No/2);
% ids=speye(d);
% idd=speye(d^N);
% h0=zeros(d,d);
% for i=1:2*j
%     m=j-i;
%     hm=(1/2)*(((j-m)*(j+m+1))^(1/2));
%     h0(i,i+1)=hm;
%     h0(i+1,i)=hm;
%     h0(d-i+1,d-i)=hm;
%     h0(d-i,d-i+1)=hm;
% end
% h=sparse(h0);
% H0=sparse(d^N,d^N);
% for i=1:N
%     if(i==1)
%         p1=h;
%     else
%         p1=ids;
%     end
%     for k=2:N
%         if(k==i)
%             p1=kron(p1,h);
%         else
%             p1=kron(p1,ids);
%         end
%     end
%     H0=H0+p1;
% end
% R=zeros(d,d);
% L=zeros(d,d);
% for i=1:2*j
%     R(i,i)=1;
%     L(d-i+1,d-i+1)=1;
% end
% r=sparse(R);
% l=sparse(L);
% Pb=idd;
% for i=1:N
%     if(i==1)
%         p4=r;
%     else
%         if(i==N)
%             p4=l;
%         else
%             p4=ids;
%         end
%     end
%     for j=2:N
%         if(j==i)
%             p4=kron(p4,r);
%         else
%             if(j==i+1)
%                 p4=kron(p4,l);
%             else
%                 p4=kron(p4,ids);
%             end
%         end
%     end
%     Pb=Pb*(idd-p4);
% end
% P=sparse(d^N,d^N);
% for i1=1:d
%     for i2=1:d
%         for i3=1:d
%             for i4=1:d
%                 for i5=1:d
%                     for i6=1:d
%                         P=P+kron(ids(:,i1),kron(ids(:,i2),kron(ids(:,i3),kron(ids(:,i4),kron(ids(:,i5),ids(:,i6))))))*kron(ids(i2,:),kron(ids(i3,:),kron(ids(i4,:),kron(ids(i5,:),kron(ids(i6,:),ids(i1,:))))));
%                     end
%                 end
%             end
%         end
%     end
% end
% Pt=(idd+P+P^2+P^3+P^4+P^5)/N;
% H=Pb*Pt*H0*Pt'*Pb';
% [Vn,En,num]=eigsub(H,Pb*Pt,800);
% save('vn12spin1','Vn');
% save('en12spin1','En');
% save('dim12spin1','num');
% st1=[1;0;0;0;0];
% E=zeros(num,1);
% ove=zeros(num,1);
% psi0=kron(st1,kron(st1,kron(st1,kron(st1,kron(st1,st1)))));
% for i=1:num
%     E(i,1)=En(i,i);
%     ove(i,1)=log10(abs(psi0'*Vn(:,i))^2);
% end
% x=-7:0.1:7;
% y=-1.2-0.089*x.^2;
% scatter(E,ove)
% hold on
% plot(x,y);
% xlim([-7,7])
% ylim([-7.5,0])

% snum=zeros(24,1);
% nsnum=zeros(num-24,1);
% k=0;
% kn=0;
% for i=1:num
%     if ove(i,1)>-1.2-0.089*E(i,1)^2
%         k=k+1;
%         snum(k,1)=i;
%     else
%         kn=kn+1;
%         nsnum(kn,1)=i;
%     end
% end

po = [];
method = [];
radius = [];
N = [];
n = [];
ms = 25;
scatplot(real(E(nsnum)),ove(nsnum),method,radius,N,n,po,ms);
hold on
scatter(real(E(snum)),ove(snum),40,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[0 0 1],'LineWidth',2)
ylim([-7.5,0])
xlim([-6.5,6.5])
xticks(-6:4:6)
yticks(-6:3:0)
% xlabel('E','FontName','Arial')
% ylabel('log_{10}|\langle E_i|\psi(0)\rangle|^2','FontName','Arial')
box on
colorbar off
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 500, 450, 400]);
% scatter(E(snum),ove(snum))
% hold on
% scatter(E(nsnum),ove(nsnum))
% hold on
% xlim([-7,7])
% ylim([-7.5,0])

% save('Snum','snum')
% save('Nsnum','nsnum')

