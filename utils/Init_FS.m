function [W1,H1,W2,H2]=Init_FS(X,ranks,noloopsflag,sparsityflag)
%% Init_FS:
% It is an initialization for the Hadamard decomposition problem, that is:
% given a matrix X of size m-by-n, compute two matrices W and H
% such that:
%   -X~=W*H'where W is m-by-(r1*r2) and H is n-by-(r1*r2),
%   -W~=face_split(W1,W2) and H~=face_split(H1,H2),
% where ranks=[r1,r2].
% The function assumes that min(m,n)>r1*r2 and it implements the
% initialization 'FS'. See HadDec_init for more details.

    r1=ranks(1);
    r2=ranks(2);
    rsvd=r1*r2;
    if noloopsflag
        proj_Bmr=@(U) projBmr_noloop(U,r1,r2);
    else
        proj_Bmr=@(U) projBmr_loop(U,r1,r2);
    end
    if sparsityflag
        [U,S,V]=svds(X,rsvd);
    else
        [U,S,V]=svd(X,'econ');
    end
    U=U(:,1:rsvd)*sqrt(S(1:rsvd,1:rsvd));
    V=V(:,1:rsvd)*sqrt(S(1:rsvd,1:rsvd));
    [W1,W2]=proj_Bmr(U);
    [H1,H2]=proj_Bmr(V);
    
end
