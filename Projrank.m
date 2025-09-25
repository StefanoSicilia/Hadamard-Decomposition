function Y=Projrank(X,W)
    
    WW=W*W';
    Y=WW*X+X*WW-WW*X*WW;

end