function [X,U,Y,V]=Had_init_flex(A)
%% Had_init_flex: 
% Given a matrix A of size nxm, it computes two matrices W and H
% such that:
%   -A~W*H'where W is nx(r^2) and H is mx(r^2),
%   -W~face_split(X,U) and H~face_split(Y,V).
% It is an inizialization for the Hadamard decomposition of A.
% The problem cannot be solved exactly if rank(A)>r^2 and this function
% also assumes that min(m,n)>r^2.

    [W,S,H]=svd(A);
    W=W(:,1:r^2)*sqrt(S(1:r^2,1:r^2));
    H=H(:,1:r^2)*sqrt(S(1:r^2,1:r^2));
    [n,m]=size(A);
    X=zeros(n,r);
    Y=zeros(m,r);
    U=zeros(n,r);
    V=zeros(m,r);
    for i=1:n
        [P,Sigma,Q]=svd(reshape(W(i,:),r,r));
        s=sqrt(Sigma(1,1));
        X(i,:)=s*P(:,1)';
        U(i,:)=s*Q(:,1)';
    end
    for i=1:m
        [F,Theta,G]=svd(reshape(H(i,:),r,r));
        t=sqrt(Theta(1,1));
        Y(i,:)=t*F(:,1)';
        V(i,:)=t*G(:,1)';
    end
end