function [W1,H1,W2,H2]=zero_rows_pert(W1,H1,W2,H2,theta)
%% zero_rows_pert: Perturbation of zero rows
% Replaces the zero rows of the matrices W1,H1,W2 and H2 with a random
% vector multiplied by theta.

    r=size(W1,2);
    rng(1)
    vW1=vecnorm(W1,2,2)==0; W1(vW1,:)=theta*rand(sum(vW1),r);
    vH1=vecnorm(H1,2,2)==0; H1(vH1,:)=theta*rand(sum(vH1),r);
    vW2=vecnorm(W2,2,2)==0; W2(vW2,:)=theta*rand(sum(vW2),r);
    vH2=vecnorm(H2,2,2)==0; H2(vH2,:)=theta*rand(sum(vH2),r);
    %disp(sum([vW1' vH1' vW2' vH2']))

end