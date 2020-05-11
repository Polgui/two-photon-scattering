function [Gamma,Lambda,Error_tracker]=Shift_site1_to_siteN(Gamma,Lambda,site1,siteN,maxSchmidtrank, Eps)
% [Gamma,Lambda,Error_tracker]=Shift_site1_to_siteN(Gamma,Lambda,site1,siteN,maxSchmidtrank, Eps)
%
% Using a series of successive swaps, brings the state of site site1 to
% site siteN. Gamma is the cell array containing the tensors for each site,
% and Lambda the array containing the singular values vector for each link.
% Error_tracker is a 2-vector whose first element is the product of the
% projection norm resulting from each swap and whose second element is the
% minimum of the projection norms.

Error_tracker=ones(1,2);
M=length(Gamma);
projection_norm=ones(abs(siteN-site1),1);
if siteN>site1
    l=1;
    for k=site1:siteN-1;
        if k==1
            [Gamma{k},Lambda{k},Gamma{k+1},projection_norm(l)]=Swap(1,Gamma{k},Lambda{k},Gamma{k+1},Lambda{k+1},maxSchmidtrank, Eps);
        elseif k==M-1
            [Gamma{k},Lambda{k},Gamma{k+1},projection_norm(l)]=Swap(Lambda{k-1},Gamma{k},Lambda{k},Gamma{k+1},1,maxSchmidtrank, Eps);
        else
            [Gamma{k},Lambda{k},Gamma{k+1},projection_norm(l)]=Swap(Lambda{k-1},Gamma{k},Lambda{k},Gamma{k+1},Lambda{k+1},maxSchmidtrank,Eps);
        end
        l=l+1;
    end
elseif site1>siteN
    swapsites=site1-1:-1:siteN;
    NumberSwaps=site1-siteN;
    for m=1:NumberSwaps;
    k=swapsites(m);    
        if k==1
            [Gamma{k},Lambda{k},Gamma{k+1},projection_norm(m)]=Swap(1,Gamma{k},Lambda{k},Gamma{k+1},Lambda{k+1},maxSchmidtrank,Eps);
        elseif k==M-1
            [Gamma{k},Lambda{k},Gamma{k+1},projection_norm(m)]=Swap(Lambda{k-1},Gamma{k},Lambda{k},Gamma{k+1},1,maxSchmidtrank,Eps);
        else
            [Gamma{k},Lambda{k},Gamma{k+1},projection_norm(m)]=Swap(Lambda{k-1},Gamma{k},Lambda{k},Gamma{k+1},Lambda{k+1},maxSchmidtrank,Eps);
        end
    end
end

if(numel(projection_norm)~=0)
    projection_norm_min=min(projection_norm(:));
    projection_norm_prod=prod(projection_norm(:));
    Error_tracker=[projection_norm_prod,projection_norm_min];
end