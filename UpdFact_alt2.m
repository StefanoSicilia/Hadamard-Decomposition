function U=UpdFact_alt2(A,Y,V,X)

    n=size(A,1);
    r=size(Y,2);
    U=zeros(n,r);
    M=(pinv(face_split(Y,V))*A')';
    for j=1:size(M,1)
        U(j,:)=pinv(X(j,:)')*reshape(M(j,:)',r,r);
    end

end