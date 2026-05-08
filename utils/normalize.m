function [A,B,alpha]=normalize(A,B,flag)
%% normalize:
% Scales the 2-norm of the columns of A and B, so that B has unit norm
% columns and A*B' is preserved.

    if flag
        alpha=vecnorm(B)';
        alpha(alpha<1e-15)=1;
        B=B*diag(1./alpha);
        A=A*diag(alpha);
    else
        alpha=ones(size(B,2),1);
    end

end