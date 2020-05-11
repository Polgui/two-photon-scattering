function rho = single_site_reduced_state_Gamma_efficient(Gamma,Lambda,site)
% rho = single_site_reduced_state_Gamma_efficient(Gamma,Lambda,site)
% returns the reduced density matrix of a given site.

    if site==1
        LGL=Lambda_multiplication(Lambda_multiplication(Gamma{site},Lambda{site},2),1,1);
        rho=tensor_contraction(LGL,conj(LGL),[1,2],[1,2]);
    elseif site==length(Gamma)
        LGL=Lambda_multiplication(Lambda_multiplication(Gamma{site},1,2),Lambda{site-1},1);
        rho=tensor_contraction(LGL,conj(LGL),[1,2],[1,2]);
    else
        LGL=Lambda_multiplication(Lambda_multiplication(Gamma{site},Lambda{site},2),Lambda{site-1},1);
        rho=tensor_contraction(LGL,conj(LGL),[1,2],[1,2]);
    end
end