% Generate a nonnegative matrix mxr with non-zero rows of density 
% alpha in (0,1) 
% Note: there is a bias towards denser matrices since we impose non-zero
% rows 

function W = sparserand(m,r,alpha)

if alpha <= 0 || alpha > 1
    error('alpha has to be in (0,1].');
end
if alpha < 1/r
    warning('Since we impose non-zero rows, the required density is too small: alpha < 1/r.');
end
% Construct non-zero pattern 
W = rand(m,r); 
W(W < 1-alpha) = 0; 
maxW = max(W'); % max per row 
findzerorow = find(maxW == 0); 
while ~isempty(findzerorow)
    W(findzerorow,:) = rand(length(findzerorow),r);
    W(W < 1-alpha) = 0; 
    maxW = max(W'); % max per row 
    findzerorow = find(maxW == 0); 
end
nnZW = find(W>0); 
W(nnZW) = rand(length(nnZW),1); 