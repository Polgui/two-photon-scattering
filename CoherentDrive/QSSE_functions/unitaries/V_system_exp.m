function UMPS = V_system_exp( gamma, imp, gammap, myphi, Omega, Delta,dimB,dt )

dimS=3;
smL=[0 0 1;
     0 0 0;
     0 0 0];
smR=[0 1 0;
     0 0 0;
     0 0 0];

a=diag(sqrt(1:dimB-1),1);   % Ladder
idB=eye(dimB);

H=KroneckerProduct(idB,idB,idB,idB,-Delta(1)*(smL'*smL)-Delta(2)*(smR'*smR)-(Omega(1)/2)*smL-((Omega(1)/2)*smL)'-(Omega(2)/2)*smR-((Omega(2)/2)*smR)');

H=H+sqrt((1-imp(1))*gamma(1)/dt)*(1i*KroneckerProduct(idB,a',idB,idB,smL)-1i*KroneckerProduct(idB,a,idB,idB,smL'))...
   +sqrt((1-imp(2))*gamma(2)/dt)*(1i*KroneckerProduct(a',idB,idB,idB,smR*exp(1i*myphi))-1i*KroneckerProduct(a,idB,idB,idB,smR'*exp(-1i*myphi)))...
   +sqrt(imp(2)*gamma(2)/dt)*(1i*KroneckerProduct(idB,a',idB,idB,smR)-1i*KroneckerProduct(idB,a,idB,idB,smR'))...
   +sqrt(imp(1)*gamma(1)/dt)*(1i*KroneckerProduct(a',idB,idB,idB,smL*exp(1i*myphi))-1i*KroneckerProduct(a,idB,idB,idB,smL'*exp(-1i*myphi)));

H=H+sqrt(gammap/dt)*(1i*KroneckerProduct(idB,idB,a',idB,smL)-1i*KroneckerProduct(idB,idB,a,idB,(smL)'));
H=H+sqrt(gammap/dt)*(1i*KroneckerProduct(idB,idB,idB,a',smR)-1i*KroneckerProduct(idB,idB,idB,a,(smR)'));
UMPS=expm(-1i.*H*dt);
UMPS=reshape(UMPS,[dimS,dimB,dimB,dimB, dimB,dimS,dimB,dimB,dimB,dimB]); 
UMPS=permute(UMPS,[(5:-1:1),5+(5:-1:1)]); % Index order: Feedback,In,System 
    

end

