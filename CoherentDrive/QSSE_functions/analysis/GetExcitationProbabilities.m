function [exc_ev]=GetExcitationProbabilities(rhoS,exc_ev)

for j=1:size(exc_ev,1)
    for k=2:size(exc_ev,2)
        exc_ev(j,k) = real(rhoS{1,k-1}(j,j));
    end
end
