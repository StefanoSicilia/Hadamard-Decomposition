function G=grad_eval(P,X)
%% grad_eval:
% Gradient of f(X_1,X_2)=0.5*norm(X-X_1.*X_2,'fro') used by 'Manopt' method
% in HadDec.

    X1=prodsvd(P.X1);
    X2=prodsvd(P.X2);
    R=X-X1.*X2;
    G=struct('X1',-R.*X2,'X2',-R.*X1);

end