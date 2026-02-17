function H2=UpdFact(X,W1,H1,W2)
%% UpdFact: update for BCD method in HadDec
% Updating function for the BCD method. For the details see  
% [WVG25] S. Wertz, A. Vandaele, N.Gillis, 
% Efficient algorithms for the Hadamard decomposition, 2025.

    H2=zeros(size(H1));
    for j=1:size(H2,1)
        H2(j,:)=hadLS(W1*H1(j,:)',W2,X(:,j))';
    end

end