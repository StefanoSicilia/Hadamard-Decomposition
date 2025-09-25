function V=UpdFact(A,X,Y,U,V)

    for j=1:size(V,1)
        V(j,:)=hadLS(X*Y(j,:)',U,A(:,j))';
    end

end