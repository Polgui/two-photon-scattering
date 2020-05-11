function [Gamma1,Lambda1,Gamma2,projection_norm]=Swap(Lambda0,Gamma1,Lambda1,Gamma2,Lambda2,maxSchmidtrank, Eps)
% [Gamma1,Lambda1,Gamma2,projection_norm]=Swap(Lambda0,Gamma1,Lambda1,Gamma2,Lambda2,maxSchmidtrank, Eps) 
%
% Swaps the states of the two adjacent sites described
% by tensors Gamma1 and Gamma2. Lambda0, Lambda1 and Lambda2 are the
% singular values vectors, respectively for the bond on the left of site 1, 
% between sites 1 and 2, and on the right of site 2.
%
% The first/second dimension of each Gamma tensor should be the length of
% the bond on the left/right of the corresponding site. The dimension of
% the third index is the physical dimension of the site.
%
% projection_norm is the 2-norm of the new Lambda1. 
% This vector is truncated to the first maxSchmidtrank highest values. 


Bondlength0=size(Gamma1,1);     % dimension of the left link
Bondlength2=size(Gamma2,2);     % dimension of the right link
dim1=size(Gamma1,3);            % physical dimension of the first site
dim2=size(Gamma2,3);            % physical dimension of the second site

% Merge all tensors into one tensor of dimension (Bondlength0, dim1, Bondlength2, dim2)
A1=Lambda_multiplication(Gamma1,Lambda0,1);
B2=Lambda_multiplication(Gamma2,Lambda2,2);
C1=Lambda_multiplication(A1,Lambda1,2);
C12=tensor_contraction(C1,B2,2,1);

% Swap by reordering the indices into (Bondlength0, dim2, Bondlength2, dim1)
C21=permute(C12,[1,4,3,2]);

%% Perform the new SVD

% Reshape C21 into a matrix and apply the SVD
C21=reshape(C21,[Bondlength0*dim2,Bondlength2*dim1]); 

[U,L1,V,projection_norm]=MySVD(C21, maxSchmidtrank, Eps);
Bondlength1=length(L1);

% Reshape the resulting U and V matrices into the original tensors
A1=permute(reshape(U,[Bondlength0,dim2,Bondlength1]),[1,3,2]);
B2=reshape(V',[Bondlength1,Bondlength2,dim1]);
Gamma1=Lambda_multiplication(A1,1./Lambda0,1);
Gamma2=Lambda_multiplication(B2,1./Lambda2,2);
Lambda1=L1;

end
    