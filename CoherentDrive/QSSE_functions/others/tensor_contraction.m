function C=tensor_contraction(A,B,ContractionsA,ContractionsB) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contracts the tensors A and B along the dimensions specified in ContractionsA and ContractionsB
%
% Simple example: let A be a NxM tensor (that is a Marix) and B a MxK
% tensor, then C=tensor_contraction(A,B,2,1) is a NxK tensor given by the
% matrix multiplication C=A*B.
%
% More elaborate example: let A be a NxMxPxQ tensor and B be a RxSxQxTxP tensor. 
% Then C=tensor_contraction(A,B,[3,4],[5,3]) gives a NxMxRxSxT tensor.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rankA=length(size(A));
rankB=length(size(B));
rankContractionsA=length(ContractionsA);
rankContractionsB=length(ContractionsB);

ContractionsA=ContractionsA(ContractionsA<=rankA);
ContractionsB=ContractionsB(ContractionsB<=rankB);


NonContractionsA=[1:rankA];
NonContractionsA(ContractionsA)=[];
NonContractionsB=[1:rankB];
NonContractionsB(ContractionsB)=[];

sizeA=size(A);
sizeContractionsA=sizeA(ContractionsA);
sizeB=size(B);
sizeContractionsB=sizeB(ContractionsB);
sizeNonContractionsA=sizeA;
sizeNonContractionsA(ContractionsA)=[];
sizeNonContractionsB=sizeB;
sizeNonContractionsB(ContractionsB)=[];

permutevecA=cat(2,NonContractionsA,ContractionsA);
AA=permute(A,permutevecA);

permutevecB=cat(2,ContractionsB,NonContractionsB);

BB=permute(B,permutevecB);

AA=reshape(AA,[prod(sizeNonContractionsA),prod(sizeContractionsA)]);

BB=reshape(BB,[prod(sizeContractionsB),prod(sizeNonContractionsB)]);

CC=AA*BB;

if and(rankContractionsA-rankA==0,rankContractionsB-rankB==0)
    sizeNonContractionsA=1;
    sizeNonContractionsB=1;
elseif and(rankContractionsA==rankA,rankContractionsB+1==rankB)
    sizeNonContractionsA=1;
elseif and(rankContractionsA+1==rankA,rankContractionsB==rankB)
    sizeNonContractionsB=1;
end

reshapevec=cat(2,sizeNonContractionsA,sizeNonContractionsB);
C=reshape(CC,reshapevec);
