function H2=UpdFact2(X,W1,H1,W2)
%% UpdFact2: update for BCD method in HadDec
% Same as UpdFact, but avoids the use of hadLS for each row of the matrix.
% For the details see [WVG25] S. Wertz, A. Vandaele, N.Gillis, 
% Efficient algorithms for the Hadamard decomposition, 2025.
% Same as UpdFact, but it uses Matlab reshapings to avoid loops.

    [m,n]=size(X);
    X1=W1*H1';

    % All-at-once operations for the rows of H2
    H=pagemtimes(permute(W2.*reshape(X1.^2,m,1,n),[2 1 3]),W2);
    D=pagemtimes(permute(W2,[2 1]),reshape(X1.*X,m,1,n)); 
    H2=pagemldivide(H,D); % least square solution
    H2=squeeze(permute(H2,[3 1 2]));

end