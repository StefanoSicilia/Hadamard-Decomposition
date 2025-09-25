function Z=column_Had(X,U)
    
    [~,r]=size(X);
    X=[1:r;X];
    X=repmat(X,1,r);
    X=sortrows(X');
    X=X(:,2:end)';
    U=repmat(U,1,r);
    Z=X.*U;

end