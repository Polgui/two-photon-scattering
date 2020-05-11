function S=EntropyProfile(Lambda)
%  S=EntropyProfile(Lambda) returns the entanglement entropy arising from a
%  single cut of the whole system, where the cut varies between each
%  component of Lambda.
M=length(Lambda);
S=zeros(1,M);
for k=1:M
    S(k)=-sum(Lambda{k}.^2.*log(Lambda{k}.^2));
end

end

