function [A,B,alpha]=normalize(A,B)
%% normalize:
% Scales the 2-norm of the columns of A and B, so that B has unit norm
% columns and A*B' is preserved.

    alpha=vecnorm(B)';
    alpha(alpha<1e-15)=1;
    B=B*diag(1./alpha);
    A=A*diag(alpha);

end