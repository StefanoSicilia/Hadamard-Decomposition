function [W1,H1,W2,H2,info]=Init_best(X,r,noloopsflag)
%% Init_best: Initialization for Hadamard decomposition
% Computes an initialization for the Hadamard Decomposition of X based on 4 
% possible choices, inspired by X=X1.*X2, X=W*H' with W and H with rows
% that are vectorization of rank-1 matrix and a rank-r^2 SVD X=U*V'.
%
% 1) SVD-based (Wertz et al. 2025) X1=sqrt(abs(X)) X2=sgn(X).*X1
% 2) Face Splitting (FS): W=proj_Bmr(U), H=proj_Bmr(V)
% 3) FS-Left: H=proj_Bmr(V), W=proj_Bmr((pinv(H)*X')')
% 4) FS-right: W=proj_Bmr(U), H=proj_Bmr(pinv(W)*X)
%
% It selects the best one according to the smallest error.
%
% Note: if r^2>min(size(A)), only SVD-based initialization is well defined
% and hence just that will be considered.

    [m,n]=size(X);

    % 1) SVD-based initialization
    M=sqrt(abs(X));
    N=sign(X).*M;
    [U1,S1,V1]=svd(M);
    W1temp=U1(:,1:r)*sqrt(S1(1:r,1:r));
    H1temp=V1(:,1:r)*sqrt(S1(1:r,1:r));
    [U2,S2,V2]=svd(N);
    W2temp=U2(:,1:r)*sqrt(S2(1:r,1:r));
    H2temp=V2(:,1:r)*sqrt(S2(1:r,1:r));
    err_SVD=norm(X-(W1temp*H1temp').*(W2temp*H2temp'),'fro');
    err=err_SVD;
    method='SVD-based';
    W1=W1temp;
    W2=W2temp;
    H1=H1temp;
    H2=H2temp;

    if noloopsflag
        proj_Bmr=@(U) projBmr_noloop(U);
    else
        proj_Bmr=@(U) projBmr_loop(U);
    end

    %% Rank-r^2 SVD
    if r^2<=min([m,n])
        [U,S,V]=svd(X);
        U=U(:,1:r^2)*sqrt(S(1:r^2,1:r^2));
        V=V(:,1:r^2)*sqrt(S(1:r^2,1:r^2));

        % Face splitting initializations common matrices
        [W1temp,W2temp]=proj_Bmr(U);
        [H1temp,H2temp]=proj_Bmr(V);
    
        % 2) FS
        err_FS=norm(X-(W1temp*H1temp').*(W2temp*H2temp'),'fro');
        if err_FS<err
            err=err_FS;
            method='FS';
            W1=W1temp;
            W2=W2temp;
            H1=H1temp;
            H2=H2temp;
        end
    
        % 3) FS-Left
        [W1temp2,W2temp2]=proj_Bmr((pinv(face_split(H1temp,H2temp))*X')');
        err_FSL=norm(X-(W1temp2*H1temp').*(W2temp2*H2temp'),'fro');
        if err_FSL<err
            err=err_FSL;
            method='FSL';
            W1=W1temp2;
            W2=W2temp2;
            H1=H1temp;
            H2=H2temp;
        end

        % 4) FS-Right
        [H1temp2,H2temp2]=proj_Bmr((pinv(face_split(W1temp,W2temp))*X)');
        err_FSR=norm(X-(W1temp*H1temp2').*(W2temp*H2temp2'),'fro');
        if err_FSR<err
            err=err_FSR;
            method='FSR';
            W1=W1temp;
            W2=W2temp;
            H1=H1temp2;
            H2=H2temp2;
        end
    end

    % Final selection
    info=struct('method',method,'err',0.5*err^2);
    
end