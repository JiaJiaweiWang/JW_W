load('vn9spin1','Vn');
load('dim9spin1','num');
j=1;
d=round(2*j+1);
N=9;
subnum=2;
Vsub=zeros(d^subnum,d^subnum,num);
for i=1:num
    Vsub(:,:,i)=ptrace(Vn(:,i)*Vn(:,i)',d,N,subnum);
end
save('vsub9spin1','Vsub')

function [Msub]=ptrace(M,d,N,n)
idd=speye(d,d);
idA=idd;
idB=idd;
for i=1:n-1
    idA=kron(idA,idd);
end
for i=1:N-n-1
    idB=kron(idB,idd);
end
Msub=zeros(d^n,d^n);
for k=1:d^(N-n)
    Msub=Msub+kron(idA,idB(k,:))*M*kron(idA,idB(:,k));
end
end