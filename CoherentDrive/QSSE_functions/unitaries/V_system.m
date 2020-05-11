function UMPS = V_system( gamma, Omega, Delta,dimB,dt )

dimS=3;
smL=[0 0 1;
     0 0 0;
     0 0 0];
smR=[0 1 0;
     0 0 0;
     0 0 0];

a=diag(sqrt(1:dimB-1),1);   % Ladder
dimB
idB=eye(dimB);

H=KroneckerProduct(idB,idB,-Delta(1)*(smL'*smL)-Delta(2)*(smR'*smR)-(Omega(1)/2)*smL-((Omega(1)/2)*smL)'-(Omega(2)/2)*smR-((Omega(2)/2)*smR)');

H=H+sqrt(gamma(1)/dt)*(1i*KroneckerProduct(idB,a',smL)-1i*KroneckerProduct(idB,a,smL'))...
   +sqrt(gamma(2)/dt)*(1i*KroneckerProduct(a',idB,smR)-1i*KroneckerProduct(a,idB,smR'));

UMPS=expm(-1i.*H*dt);
UMPS=reshape(UMPS,[dimS,dimB,dimB,dimS,dimB,dimB]); 
UMPS=permute(UMPS,[(3:-1:1),3+(3:-1:1)]); % Index order: Feedback,In,System 
    

end

