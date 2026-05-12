function [Z1,Z2]=UpdManifold_noloop(Z1,Z2,G,h)
%% UpdManifold_noloop: update for manBCD method in HadDec
% Updates the pair (Z1,Z2), which represents either (W1,W2) or (H1,H2).
% It follows the gradient system showed in the paper.
% Same as UpdManifold_loop, but it uses Matlab reshapings to avoid the loop.


    r1=size(Z1,2);
    r2=size(Z2,2);
    rho1=vecnorm(Z1,2,2);
    rho2=vecnorm(Z2,2,2);
    rho=rho1.*rho2;
    rho1(rho1<1e-15)=1;
    rho2(rho2<1e-15)=1;
    X=Z1./rho1;
    Y=Z2./rho2;
    idx=rho>1e-15;
    if any(idx)
        G3=reshape(G(idx,:)',r2,r1,[]);
        Xp=permute(X(idx,:),[2 3 1]);
        Yp=permute(Y(idx,:),[2 3 1]);
        Gixi=pagemtimes(G3,Xp);
        yiGixi=pagemtimes(permute(Yp,[2 1 3]),Gixi);
        GitY=pagemtimes(permute(G3,[2 1 3]),Yp); 
        sta=rho(idx)./reshape(yiGixi,[],1);
        h=min(h,0.95*sta).*(sta>0)+h.*(sta<=0);
        updater=sqrt(1-h./sta);
        rho1(idx)=rho1(idx).*updater;
        rho2(idx)=rho2(idx).*updater;
        scale=reshape(h./(rho1(idx).*rho2(idx)),1,1,[]);
        X(idx,:)=permute(Xp+scale.*(-GitY+Xp.*yiGixi),[3 1 2]);
        Y(idx,:)=permute(Yp+scale.*(-Gixi+Yp.*yiGixi),[3 1 2]);
    end
    Z1=X.*rho1;
    Z2=Y.*rho2;

end