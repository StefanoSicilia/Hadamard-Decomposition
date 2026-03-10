function [Z1,Z2]=UpdManifold_loop(Z1,Z2,G,h)
%% UpdManifold_loop: update for manBCD method in HadDec
% Updates the pair (Z1,Z2), which represents either (W1,W2) or (H1,H2).
% It follows the gradient system showed in the paper.
% Same as UpdManifold_noloop, but it uses a loop.

    [m,r]=size(Z1);
    rho1=vecnorm(Z1,2,2);
    rho2=vecnorm(Z2,2,2);
    rho=rho1.*rho2;
    rho1(rho1<1e-15)=1;
    rho2(rho2<1e-15)=1;
    X=Z1./rho1;
    Y=Z2./rho2;
    v=1:m;
    for i=v(rho>1e-15)
        Gi=reshape(G(i,:),r,r);
        Gixi=Gi*X(i,:)';
        yiGixi=Y(i,:)*Gixi;
        X(i,:)=X(i,:)+h/rho(i)*(-Gi'*Y(i,:)'+0.5*(yiGixi)*X(i,:)')';
        Y(i,:)=Y(i,:)+h/rho(i)*(-Gixi+0.5*(yiGixi)*Y(i,:)')';
    end
    Z1=X.*rho1;
    Z2=Y.*rho2;

end