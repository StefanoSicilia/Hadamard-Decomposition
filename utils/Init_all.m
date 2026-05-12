function output=Init_all(X,ranks,noloopsflag,sparsityflag)
%% Init_all: All initializations for Hadamard decomposition
% Given ranks=[r1,r2], it computes all initializations for the 
% (r1,r2)-Hadamard Decomposition of X based on 4 possible choices, 
% inspired by X=X1.*X2, X=W*H' with W and H with rows that are 
% vectorizations of rank-1 matrix and by the rank-(r1*r2) SVD X=U*V'.
%
% 1) SVD-based (Wertz et al. 2025) X1=sqrt(abs(X)) X2=sgn(X).*X1
% 2) Face Splitting (FS): W=proj_Bmr(U), H=proj_Bmr(V)
% 3) FS-Left: H=proj_Bmr(V), W=proj_Bmr((pinv(H)*X')')
% 4) FS-right: W=proj_Bmr(U), H=proj_Bmr(pinv(W)*X)
%
% Note: it works only if r1*r2>=min(size(A)).

    output=[];
    r1=ranks(1);
    r2=ranks(2);

    % 1) SVD-based initialization
    M=sqrt(abs(X));
    N=sign(X).*M;
    if sparsityflag
        [U1,S1,V1]=svds(M,r1);
        [U2,S2,V2]=svds(N,r2);
    else
        [U1,S1,V1]=svd(M,'econ');
        [U2,S2,V2]=svd(N,'econ');
    end
    W1temp=U1(:,1:r1)*sqrt(S1(1:r1,1:r1));
    H1temp=V1(:,1:r1)*sqrt(S1(1:r1,1:r1));
    W2temp=U2(:,1:r2)*sqrt(S2(1:r2,1:r2));
    H2temp=V2(:,1:r2)*sqrt(S2(1:r2,1:r2));
    output.W1_svd=W1temp;
    output.W2_svd=W2temp;
    output.H1_svd=H1temp;
    output.H2_svd=H2temp;

    if noloopsflag
        proj_Bmr=@(U) projBmr_noloop(U,r1,r2);
    else
        proj_Bmr=@(U) projBmr_loop(U,r1,r2);
    end

    rsvd=r1*r2;
    if sparsityflag
        [U,S,V]=svds(X,rsvd);
        U=U*sqrt(S);
        V=V*sqrt(S);
    else
        [U,S,V]=svd(X,'econ');
        U=U(:,1:rsvd)*sqrt(S(1:rsvd,1:rsvd));
        V=V(:,1:rsvd)*sqrt(S(1:rsvd,1:rsvd));
    end
    

    % Face splitting initializations common matrices
    [W1temp,W2temp]=proj_Bmr(U);
    [H1temp,H2temp]=proj_Bmr(V);

    % 2) FS
    output.W1_FS=W1temp;
    output.W2_FS=W2temp;
    output.H1_FS=H1temp;
    output.H2_FS=H2temp;

    % 3) FS-Left
    [W1temp2,W2temp2]=proj_Bmr((pinv(face_split(H1temp,H2temp))*X')');
    output.W1_FSL=W1temp2;
    output.W2_FSL=W2temp2;
    output.H1_FSL=H1temp;
    output.H2_FSL=H2temp;

    % 4) FS-Right
    [H1temp2,H2temp2]=proj_Bmr((pinv(face_split(W1temp,W2temp))*X)');
    output.W1_FSR=W1temp;
    output.W2_FSR=W2temp;
    output.H1_FSR=H1temp2;
    output.H2_FSR=H2temp2;
   
end