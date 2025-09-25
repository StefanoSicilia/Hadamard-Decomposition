function X=UpdFact_alt(A,Y,V,U)

    n=size(A,1);
    r=size(Y,2);
    X=zeros(n,r);
    M=(pinv(face_split(Y,V))*A')';
    for j=1:size(M,1)
        X(j,:)=pinv(U(j,:)')*reshape(M(j,:),r,r);
    end

end