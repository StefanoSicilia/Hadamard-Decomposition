function [W1,H1,W2,H2]=Init_SVDbased(X,r)
%% Init_SVDbased: 
% It is an initialization for the Hadamard decomposition problem, that is:
% given a matrix X of size m-by-n, compute two matrices X1 and X2
% such that X~X1*X2', where X1 and X2 are rank-r mby-n matrices.
% The function implements the initialization 'SVD-based'.

    M=sqrt(abs(X));
    N=sign(X).*M;
    [U1,S1,V1]=svd(M);
    W1=U1(:,1:r)*sqrt(S1(1:r,1:r));
    H1=V1(:,1:r)*sqrt(S1(1:r,1:r));
    [U2,S2,V2]=svd(N);
    W2=U2(:,1:r)*sqrt(S2(1:r,1:r));
    H2=V2(:,1:r)*sqrt(S2(1:r,1:r));

end