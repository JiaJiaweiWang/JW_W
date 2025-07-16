s1=sparse([0,1/(2^(1/2)),0;1/(2^(1/2)),0,0;0,0,0]);
s2=sparse([0,0,0;0,0,1/(2^(1/2));0,1/(2^(1/2)),0]);
l=sparse([0,0,0;0,0,0;0,0,1]);
r=sparse([1,0,0;0,0,0;0,0,0]);
sx=s1+s2;
id3=speye(3);
d=9; %system size
idd=speye(3^d);
H=sparse(3^d,3^d);
for i=1:d
    if(i==1)
        p1=sx;
    else
        p1=id3;
    end
    for j=2:d
        if(j==i)
            p1=kron(p1,sx);
        else
            p1=kron(p1,id3);
        end
    end
    H=H+p1;
end
for i=1:d
    if(i==1)
        p2=s1;
    else
        if(i==d)
            p2=l;
        else
            p2=id3;
        end
    end
    for j=2:d
        if(j==i)
            p2=kron(p2,s1);
        else
            if(j==i+1)
                p2=kron(p2,l);
            else
                p2=kron(p2,id3);
            end
        end
    end
    H=H-p2;
end
for i=1:d
    if(i==1)
        p3=r;
    else
        if(i==d)
            p3=s2;
        else
            p3=id3;
        end
    end
    for j=2:d
        if(j==i)
            p3=kron(p3,r);
        else
            if(j==i+1)
                p3=kron(p3,s2);
            else
                p3=kron(p3,id3);
            end
        end
    end
    H=H-p3;
end
P=sparse(3^d,3^d);
for i1=1:3
    for i2=1:3
        for i3=1:3
            for i4=1:3
                for i5=1:3
                    for i6=1:3
                        for i7=1:3
                            for i8=1:3
                                for i9=1:3
                                    P=P+kron(id3(:,i1),kron(id3(:,i2),kron(id3(:,i3),kron(id3(:,i4),kron(id3(:,i5),kron(id3(:,i6),kron(id3(:,i7),kron(id3(:,i8),id3(:,i9)))))))))*kron(id3(i2,:),kron(id3(i3,:),kron(id3(i4,:),kron(id3(i5,:),kron(id3(i6,:),kron(id3(i7,:),kron(id3(i8,:),kron(id3(i9,:),id3(i1,:)))))))));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
Pt=(idd+P+P^2+P^3+P^4+P^5+P^6+P^7+P^8)/d;
Po=idd;
for i=1:d
    if(i==1)
        p4=r;
    else
        if(i==d)
            p4=l;
        else
            p4=id3;
        end
    end
    for j=2:d
        if(j==i)
            p4=kron(p4,r);
        else
            if(j==i+1)
                p4=kron(p4,l);
            else
                p4=kron(p4,id3);
            end
        end
    end
    Po=Po*(idd-p4);
end
[Vn,En,num]=eigsub(H,Pt*Po,1750);
save('vn9spin1','Vn');
save('en9spin1','En');
save('dim9spin1','num');