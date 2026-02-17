function [Z1,Z2]=Upd_manifold(Z1,Z2,G,h)
%% Upd_manifold: update for manBCD method in HadDec
% Updates the pair (Z1,Z2), which represents either (W1,W2) or (H1,H2).
% It follows the gradient system showed in the paper.

    [m,r]=size(Z1);
    
    for i=1:m
        Gi=reshape(G(i,:),r,r);
        Gixi=Gi*Z1(i,:)';
        yiGixi=Z2(i,:)*Gixi;
        Z1(i,:)=Z1(i,:)+h*(-Gi'*Z2(i,:)'+0.5*(yiGixi)*Z1(i,:)')';
        Z2(i,:)=Z2(i,:)+h*(-Gixi+0.5*(yiGixi)*Z2(i,:)')';
    end

end