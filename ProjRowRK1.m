function A=ProjRowRK1(A)

    n=size(A,1);
    r=sqrt(n);
    for i=1:n
        [P,Sigma,Q]=svd(reshape(A(i,:),r,r));
        A(i,:)=reshape(Sigma(1,1)*P(:,1)*Q(:,1)',r^2,1);
    end

end