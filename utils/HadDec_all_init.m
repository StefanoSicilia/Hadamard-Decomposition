function output=HadDec_all_init(X,r)
%% HadDec_all_init: All initializations for Hadamard decomposition
% Computes all initializations for the Hadamard Decomposition of X based on 
% 4 possible choices, inspired by X=X1.*X2, X=W*H' with W and H with rows
% that are vectorization of rank-1 matrix and a rank-r^2 SVD X=U*V'.
%
% 1) SVD-based (Wertz et al. 2025) X1=sqrt(abs(X)) X2=sgn(X).*X1
% 2) Face Splitting (FS): W=proj_Bmr(U), H=proj_Bmr(V)
% 3) FS-Left: H=proj_Bmr(V), W=proj_Bmr((pinv(H)*X')')
% 4) FS-right: W=proj_Bmr(U), H=proj_Bmr(pinv(W)*X)
%
% Note: it works only if r^2>min(size(A)).

    output=[];

    % 1) SVD-based initialization
    M=sqrt(abs(X));
    N=sign(X).*M;
    [U1,S1,V1]=svd(M);
    W1temp=U1(:,1:r)*sqrt(S1(1:r,1:r));
    H1temp=V1(:,1:r)*sqrt(S1(1:r,1:r));
    [U2,S2,V2]=svd(N);
    W2temp=U2(:,1:r)*sqrt(S2(1:r,1:r));
    H2temp=V2(:,1:r)*sqrt(S2(1:r,1:r));
    output.W1_svd=W1temp;
    output.W2_svd=W2temp;
    output.H1_svd=H1temp;
    output.H2_svd=H2temp;

    [U,S,V]=svd(X);
    U=U(:,1:r^2)*sqrt(S(1:r^2,1:r^2));
    V=V(:,1:r^2)*sqrt(S(1:r^2,1:r^2));

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