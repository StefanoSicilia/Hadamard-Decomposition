function [W1,H1,W2,H2]=Init_FSL(X,r)
%% Init_FSL: 
% It is an initialization for the Hadamard decomposition problem, that is:
% given a matrix X of size m-by-n, compute two matrices W and H
% such that:
%   -X~W*H'where W is m-by-(r^2) and H is n-by-(r^2),
%   -W~face_split(W1,W2) and H~face_split(H1,H2).
% The function assumes that min(m,n)>r^2 and it implements the
% initialization 'FSL'.

    [~,S,V]=svd(X);
    V=V(:,1:r^2)*sqrt(S(1:r^2,1:r^2));
    [H1,H2]=proj_Bmr(V);
    [W1,W2]=proj_Bmr((pinv(face_split(H1,H2))*X')');

end