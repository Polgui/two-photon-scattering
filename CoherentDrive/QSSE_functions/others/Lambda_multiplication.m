function GL=Lambda_multiplication(Gamma,Lambda,index)
% Lambda_multiplcation(Gamma,Lambda,index) multiplies Gamma by Lambda
% elementwise, where Lambda is a column vector. Gamma has two MPS indices
% (left and right) and one physical index. index=1 performs the
% multiplication on the left index of Gamma, index=2 on the right one.

dims=size(Gamma);
dims(index)=[];
permutations=[2:index,1,index+1:length(size(Gamma))];
GL=Gamma.*permute(repmat(Lambda,[1,dims]),permutations);
end