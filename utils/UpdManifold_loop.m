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
        sta=rho(i)/yiGixi;
        h=min(h,0.95*sta).*(sta>0)+h.*(sta<=0);
        updater=sqrt(1-h/sta);
        rho1(i)=rho1(i).*updater;
        rho2(i)=rho2(i).*updater;
        scale=h/(rho1(i)*rho2(i));
        X(i,:)=X(i,:)+scale*(-Gi'*Y(i,:)'+yiGixi*X(i,:)')';
        Y(i,:)=Y(i,:)+scale*(-Gixi+yiGixi*Y(i,:)')';
    end
    Z1=X.*rho1;
    Z2=Y.*rho2;

end