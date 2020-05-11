function U_loss=U_WaveguideLosses(dimB,theta)
% U_loss=U_WaveguideLosses(dimB,theta) 
% Returns the unitary tensor of dimension (dimB dimB dimB dimB)
% corresponding to the interaction between a waveguide time-bin and a loss
% channel time-bin. The loss is implemented as a beam-splitter coupling the
% two bins, with theta the beam-splitter angle. theta=0 corresponds to no
% loss, theta=pi/2 to 100% losses

a=diag(sqrt(1:dimB-1),1);
H=KroneckerProduct(a',a)+KroneckerProduct(a,a');
U_loss=reshape(expm(1i*theta*H),[dimB,dimB,dimB,dimB]);

end

