function [Gamma,Lambda,ErrorTracker] = LocalUnitaryEvolution(Gamma,Lambda,siteS,UMPS,maxSchmidtrank, Eps)
% [Gamma,Lambda, ErrorTracker] = LocalUnitaryEvolution(Gamma,Lambda,siteS,UMPS,maxSchmidtrank, Eps)
%
% Performs the local unitary evolution corresponding to UMPS. Gamma and
% Lambda are the full cell arrays, and siteS the index for the system component. The
% system tensor should be located on the right of the other ones. 
% The dimension of UMPS is thus [dimB ... dimB dimS dimB ... dimB dimS] where 
% dimS is the physical dimension of the system component and dimB of the 
% photonic field.

% Physical dimensions
dimB=size(Gamma{siteS-1},3);
dimS=size(Gamma{siteS},3);

N_B=length(size(UMPS))/2 -1; %Number of time bins interacting

%% Merge everything

siteFeedback=siteS-N_B;
Gamma_tot=Gamma{siteFeedback};
for l=1:N_B
Gamma_tot=tensor_contraction(Lambda_multiplication(Gamma_tot,Lambda{siteFeedback+(l-1)},l+1)...
                                    ,Gamma{siteFeedback+l},l+1,1);
end

Gamma_tot=permute(Gamma_tot,[1 N_B+2 2:N_B+1 N_B+3]);

BondlengthL=size(Gamma_tot,1); % left MPS dimension of Gamma_feedback 
BondlengthR=size(Gamma_tot,2); % right MPS dimension of Gamma_feedback 

if siteS-N_B-1 > 0
    Gamma_tot=Lambda_multiplication(Gamma_tot,Lambda{siteS-N_B-1},1);
end

Gamma_tot=Lambda_multiplication(Gamma_tot,Lambda{siteS},2);


%% Apply the unitary

Gamma_tot=tensor_contraction(Gamma_tot,UMPS,(3:length(size(Gamma_tot))),((length(size(UMPS))/2 +1):length(size(UMPS))));
Gamma_tot=permute(Gamma_tot,[1 (3:length(size(Gamma_tot))) 2]);

%% Bring everything back in the original form
ErrorTracker=ones(1,2);
for l=1:N_B
    Gamma_tot=reshape(Gamma_tot,[BondlengthL*dimB,dimB^(N_B-l)*dimS*BondlengthR]); % reshape into a matrix

    % Apply the SVD
    [U,L1,V, projection_norm]=MySVD(Gamma_tot,maxSchmidtrank,Eps);
    ErrorTracker(1)=ErrorTracker(1)*projection_norm;
    ErrorTracker(2)=min(ErrorTracker(2),projection_norm);
    
    Bondlength=length(L1);
    Lambda{siteS-N_B+(l-1)}=L1;
    % New state
    AL_feedback=permute(reshape(U,[BondlengthL,dimB,Bondlength]),[1,3,2]);
    
    if siteS-N_B+(l-2)>0
        Gamma{siteS-N_B+(l-1)}=Lambda_multiplication(AL_feedback,1./Lambda{siteS-N_B+(l-2)},1);
    end
    
    % All the rest
    V=V(:,1:Bondlength);
    
    Gamma_tot=diag(L1)*V';
    BondlengthL=Bondlength;

end

% New state for the system
Gamma{siteS}=Lambda_multiplication(permute(reshape(V',[BondlengthL,dimS,BondlengthR]),[1,3,2]),1./Lambda{siteS},2);

