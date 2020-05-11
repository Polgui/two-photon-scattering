function UMPS = V_system( gamma, gammap,myphi,Omega, Delta,dimB,dt )

dimS=3;
smL=[0 0 1;
     0 0 0;
     0 0 0];
smR=[0 1 0;
     0 0 0;
     0 0 0];

a=diag(sqrt(1:dimB-1),1);   % Ladder
idB=eye(dimB);

H=KroneckerProduct(idB,idB,idB,-Delta(1)*(smL'*smL)-Delta(2)*(smR'*smR)-(Omega(1)/2)*smL-((Omega(1)/2)*smL)'-(Omega(2)/2)*smR-((Omega(2)/2)*smR)');

H=H+sqrt(gamma(1)/dt)*(1i*KroneckerProduct(idB,a',idB,smL)-1i*KroneckerProduct(idB,a,idB,smL'))...
   +sqrt(gamma(2)/dt)*(1i*KroneckerProduct(a',idB,idB,smR*exp(1i*myphi))-1i*KroneckerProduct(a,idB,idB,smR'*exp(-1i*myphi)))...
   +sqrt(gammap/dt)*(1i*KroneckerProduct(idB,idB,a',smR+smL)-1i*KroneckerProduct(idB,idB,a,smR'+smL'));

UMPS=expm(-1i.*H*dt);
UMPS=reshape(UMPS,[dimS,dimB,dimB,dimB,dimS,dimB,dimB,dimB]); 
UMPS=permute(UMPS,[(4:-1:1),4+(4:-1:1)]); % Index order: Feedback,In,Drive,System 
    

end

