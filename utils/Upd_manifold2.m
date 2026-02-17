function [Z1,Z2]=Upd_manifold2(Z1,Z2,G,h)
%% Upd_manifold2: update for manBCD method in HadDec
% Updates the pair (Z1,Z2), which represents either (W1,W2) or (H1,H2).
% It follows the gradient system showed in the paper.
% Same as Upd_manifold, but it avoids the loop and uses Matlab reshapings.

    [m,r]=size(Z1);
    G3=reshape(G',r,r,m);
    Z1p=permute(Z1,[2 3 1]);
    Z2p=permute(Z2,[2 3 1]);
    Gixi=pagemtimes(G3,Z1p);
    yiGixi=pagemtimes(permute(Z2p,[2 1 3]),Gixi);
    GitZ2=pagemtimes(permute(G3,[2 1 3]),Z2p);
    Z1=permute(Z1p+h*(-GitZ2+0.5*Z1p.*yiGixi),[3 1 2]);
    Z2=permute(Z2p+h*(-Gixi+0.5*Z2p.*yiGixi),[3 1 2]);

end