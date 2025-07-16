function [Vn,En,dim]=eigsub(H,Sym,d)% 基于估计的空间大小求解对称性子空间内的本征值与本征态
[Symbase,SymE]=eigs(Sym,d);
dim=d;
for i=1:d
    if SymE(i,i)<0.5
        dim=i-1;
        break
    end
end
if(dim<d)
    SymB=Symbase(:,1:dim);
    Hsub=zeros(dim,dim);
    for i=1:dim
        for k=1:dim
            Hsub(i,k)=SymB(:,i)'*H*SymB(:,k);
        end
    end
    Hsubs=sparse(Hsub);
    [Vs,En]=eigs(Hsubs,dim);
    Vn=SymB*Vs;
else
    Vn=0;
    En=0;
end
end