% Generate low-rank sparse nonnegative matrix whose NMF with density beta 

function [W,H,X] = sparserandNMF(m,n,r,beta)

if beta <= 0 || beta > 1
    error('beta has to be in (0,1].');
end
gamma = exp(log(1-beta)/r); 
alpha = sqrt(1-gamma); 
W = sparserand(m,r,alpha); 
H = sparserand(r,n,alpha); 
X = W*H; 