function H2=UpdFact_noloop(X,W1,H1,W2)
%% UpdFact_noloop: update for BCD method in HadDec
% For the details see [WVG25] S. Wertz, A. Vandaele, N.Gillis, 
% Efficient algorithms for the Hadamard decomposition, 2025.
% Same as UpdFact_loop, but it uses Matlab reshapings to avoid the loop.

    [m,n]=size(X);
    X1=W1*H1';

    % All-at-once operations for the rows of H2
    H=pagemtimes(permute(W2.*reshape(X1.^2,m,1,n),[2 1 3]),W2);
    D=pagemtimes(permute(W2,[2 1]),reshape(X1.*X,m,1,n)); 
    H2=pagemldivide(H,D); % least square solution
    H2=squeeze(permute(H2,[3 1 2]));

end