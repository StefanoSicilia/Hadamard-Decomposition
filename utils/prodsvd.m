function X=prodsvd(A)
%% prodsvd:
% Computes the explicit matrix written in the form of the rank-1 manifold
% of Manopt.

    X=A.U*A.S*A.V';

end