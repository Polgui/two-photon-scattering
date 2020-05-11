function Correlation=AutoCorrelation2(Gamma,Lambda,operator11,operator12,operator21,operator22,site1,site2)

% site 1 is the site that is fixed, site 2 is varied

operator1=operator11*operator12;
operator2=operator21*operator22;

Correlation=zeros(1,abs(site2-site1));
contr_index=(site2<site1)+1;
open_index=(site2>=site1)+1;

if site1==1
    C1=Lambda_multiplication(Gamma{site1},Lambda{site1},2);
elseif site1==length(Gamma)
    C1=Lambda_multiplication(Gamma{site1},Lambda{site1-1},1);
else
    C1=Lambda_multiplication(Lambda_multiplication(Gamma{site1},Lambda{site1},2),Lambda{site1-1},1);
end

rho=single_site_reduced_state_Gamma_efficient(Gamma,Lambda,site1);
occ1=trace(operator1*rho);

OC12=tensor_contraction(C1,operator11*operator21*operator22*operator12,3,2);
Correlation(1)=tensor_contraction(OC12,conj(C1),[1,2,3],[1,2,3])/occ1^2;

OC1=tensor_contraction(C1,operator1,3,2);
OT1=tensor_contraction(OC1,conj(C1),[contr_index,3],[contr_index,3]);

sites=(site1+sign(site2-site1):sign(site2-site1):site2);

for d=1:length(sites)
    B2=Gamma{sites(d)};
    
    rho=single_site_reduced_state_Gamma_efficient(Gamma,Lambda,sites(d));
    occ2=trace(operator2*rho);
    
    if sites(d)-(site2<site1)~=length(Gamma) && sites(d)-(site2<site1)~=0
        B2=Lambda_multiplication(B2,Lambda{sites(d)-(site2<site1)},open_index);
    end
    
    OB2=tensor_contraction(B2,operator2,3,2);
    OT2=tensor_contraction(OB2,conj(B2),[open_index,3],[open_index,3]);
    Correlation(d+1)=tensor_contraction(OT1,OT2,[1,2],[1,2])/(occ1*occ2);

    BOT1=tensor_contraction(B2,OT1,contr_index,1);
    OT1=tensor_contraction(BOT1,conj(B2),[3,2],[contr_index,3]);
end
