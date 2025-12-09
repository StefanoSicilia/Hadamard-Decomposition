function H=hess_eval(P,A,X,manifold)
%% hess_eval:
% Hessian of f(X_1,X_2)=0.5*norm(X-X_1.*X_2,'fro') used by HadDec_Manopt.
    
    X1=prodsvd(P.X1);
    X2=prodsvd(P.X2);
    Y=2*X1.*X2-X;
    B=manifold.tangent2ambient(P,A);
    dA1=prodsvd(B.X1);
    dA2=prodsvd(B.X2);
    H=struct('X1',dA1.*(X2.^2)+Y.*dA2,'X2',dA2.*(X1.^2)+Y.*dA1);

end