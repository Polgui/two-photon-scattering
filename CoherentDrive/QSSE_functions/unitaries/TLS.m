function UMPS = TLS( gamma, Omega, Delta,dimB,dt )
% UMPS = TLS( gamma, Omega, Delta,dimB,dt )
% Returns the Unitary tensor of dimension (dimB ... dimB 2 dimB ... dimB 2)
% corresponding to the interaction between a driven TLE and multiple waveguides 
% for a time dt. gamma is a vector containing the coupling with each waveguide, 
% Omega the Rabi frequency of the driving and Delta its detuning.

dimS=2;
sm=[0,1;0,0];               % Pauli sigma-
a=diag(sqrt(1:dimB-1),1);   % Ladder
idB=eye(dimB);
N_wav=length(gamma);

H=KroneckerProduct(KroneckerProduct(idB,N_wav),-Delta*sm'*sm-(Omega/2)*sm-(Omega/2)'*sm');

for j=1:N_wav
    H=H+sqrt(gamma(j)/dt)*(1i*KroneckerProduct(KroneckerProduct(idB,j-1),a',KroneckerProduct(idB,N_wav-j),sm)...
        -1i*KroneckerProduct(KroneckerProduct(idB,j-1),a,KroneckerProduct(idB,N_wav-j),sm'));
end

UMPS=expm(-1i.*H*dt);
UMPS=reshape(UMPS,[dimS,dimB*ones(1,N_wav),dimS,dimB*ones(1,N_wav)]); 
UMPS=permute(UMPS,[(N_wav+1:-1:1),N_wav+1+(N_wav+1:-1:1)]); % Index order: Feedback,In,System 
    

end

