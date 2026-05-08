function [W1,H1,W2,H2]=Init_SVDbased(X,ranks,sparsityflag)
%% Init_SVDbased: 
% It is an initialization for the Hadamard decomposition problem, that is:
% given a matrix X of size m-by-n, compute two matrices X1 and X2
% such that X~=X1*X2', where X1 and X2 are rank-r1 and rank-r2 
% m-by-n matrices (where ranks=[r1,r2]).
% The function implements the initialization 'SVD-based'.
% See HadDec_init for more details.

    r1=ranks(1);
    r2=ranks(2);

    M=sqrt(abs(X));
    N=sign(X).*M;
    if sparsityflag
        [U1,S1,V1]=svds(M,r1);
        [U2,S2,V2]=svds(N,r2);
    else
        [U1,S1,V1]=svd(M);
        [U2,S2,V2]=svd(N);
    end
    W1=U1(:,1:r1)*sqrt(S1(1:r1,1:r1));
    H1=V1(:,1:r1)*sqrt(S1(1:r1,1:r1));
    W2=U2(:,1:r2)*sqrt(S2(1:r2,1:r2));
    H2=V2(:,1:r2)*sqrt(S2(1:r2,1:r2));

end