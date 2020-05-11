function [U,L1,V, projection_norm]=MySVD(C21, maxSchmidtrank, Eps)
% [U,L1,V, projection_norm]=MySVD(C21, maxSchmidtrank, Eps)
% Performs the SVD of the matrix C21, but retains only the maxSchmidtrank
% higher singular values larger than Eps

[U,S,V]=svd(C21,'econ');
L1=diag(S);                     % Singular values vector

% Truncation of the singular values
Bondlength1=min(length(L1(L1>Eps)),maxSchmidtrank);
if Bondlength1==maxSchmidtrank
    warning('mytest:maxrank','Bond dimension saturated for MAXSCHMIDTRANK=%i. This might yield unprecise results.',maxSchmidtrank)
    warning('off','mytest:maxrank')
end
L1=L1(1:Bondlength1);
projection_norm=norm(L1);
L1=L1./norm(L1);
U=U(:,1:Bondlength1);
V=V(:,1:Bondlength1);


end

