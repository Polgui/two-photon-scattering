%% SETTINGS
clearvars
addpath('QSSE_functions/analysis/')
addpath('QSSE_functions/evolution/')
addpath('QSSE_functions/others/')
addpath('QSSE_functions/unitaries/')

gammapVec=[5e-1 1e-1 5e-2];
SpectrumResult=cell(1,length(gammapVec));
g2Function=cell(1,length(gammapVec));

showprogress=0;
showfigures=0;
saveresults=1;

for itergammaP=1:length(gammapVec)
   
gammap=gammapVec(itergammaP);

OmegaL=0;


dt=0.05;          % Timestep

TotalTime=120;        % Total simulation time
Nt=round(TotalTime/dt);    % Total number of timesteps

tau = 1;
M=round(tau/dt);

N_wav=2;                                % Number of waveguides

e_init=[1 0 0];

dimB=3;                                 % Dimension for bosonic modes = 3 (no photon, 1 photon, 2 photons) for each waveguide
dimS=3;                                 % Dimension for each atom

maxSchmidtrank=120; % Maximum Schmidt rank for the SVD 
Eps = 10^(-12);      % Singular values truncation   

% Decay rates in the waveguide

gamma=1.;
gammaL=gamma; 
gammaR=gamma;


myphi=0.1*2*pi;

OmegaR=OmegaL;

DeltaL=-0;
DeltaR=0;
    

    % Definition of Unitary

    UMPS=V_system([gammaL gammaR], gammap, myphi, [OmegaL OmegaR], [DeltaL DeltaR], dimB,dt);

    g=(1/sqrt(Nt))*ones(1,Nt);
    
    Gammavac=reshape([1 zeros(1,dimB-1)],[1,1,dimB]);
    
    GammaS=reshape(e_init,[1,1,dimS]);

    Gamma=cell(1,N_wav*M);   % Initial Gamma
    Lambda=cell(1,N_wav*M);  % Initial Lambda

    for n=1:N_wav*M
        Gamma{n}=Gammavac;
        Lambda{n}=1;
    end

    Gamma{N_wav*M+1}=GammaS;
    Lambda{N_wav*M+1}=1;
    
    siteS=N_wav*M+1;

    for j=1:Nt
        Gamma{N_wav*M+1+j}=Gammavac;
        Lambda{N_wav*M+1+j}=1;
    end
    
    for j=1:Nt-1
        Lambda{N_wav*M+1+Nt+j}=[g(j+1:end)*g(j+1:end)'; sqrt(2*g(j+1:end)*g(j+1:end)'*g(1:j)*g(1:j)'); g(1:j)*g(1:j)'];
    end
    Lambda{N_wav*M+1+Nt+Nt}=1;
    
    Gamma{N_wav*M+1+Nt+1}(:,:,1)=[1 0 0];
    Gamma{N_wav*M+1+Nt+1}(:,:,2)=[0 1 0];
    Gamma{N_wav*M+1+Nt+1}(:,:,3)=[0 0 1];
    for j=2:Nt-1
        x2=g(j)/abs(g(j)) * sqrt(1/Lambda{N_wav*M+1+Nt+j}(2)^2 - abs(g(j)^2)^2/Lambda{N_wav*M+1+Nt+j}(2)^2/Lambda{N_wav*M+1+Nt+j-1}(1)^2 ...
        -Lambda{N_wav*M+1+Nt+j}(1)^2/Lambda{N_wav*M+1+Nt+j}(2)^2/Lambda{N_wav*M+1+Nt+j-1}(1)^2);
        
        x3=g(j)/abs(g(j)) * sqrt(1/Lambda{N_wav*M+1+Nt+j-1}(2)^2 - abs(g(j)^2)^2/Lambda{N_wav*M+1+Nt+j}(3)^2/Lambda{N_wav*M+1+Nt+j-1}(2)^2 ...
        -Lambda{N_wav*M+1+Nt+j-1}(3)^2/Lambda{N_wav*M+1+Nt+j-1}(2)^2/Lambda{N_wav*M+1+Nt+j}(3)^2);
        
        Gamma{N_wav*M+1+Nt+j}(:,:,1)=[1/Lambda{N_wav*M+1+Nt+j-1}(1) 0 0; 0 sqrt(1/Lambda{N_wav*M+1+Nt+j-1}(2)^2 - ...
                                                            abs(x2)^2*Lambda{N_wav*M+1+Nt+j-1}(1)^2/Lambda{N_wav*M+1+Nt+j-1}(2)^2) 0;...
                                                            0 0 1/Lambda{N_wav*M+1+Nt+j}(3)];
        Gamma{N_wav*M+1+Nt+j}(:,:,2)=[0 x2 0; 0 0 x3; 0 0 0];
        Gamma{N_wav*M+1+Nt+j}(:,:,3)=[0 0 g(j)^2/Lambda{N_wav*M+1+Nt+j-1}(1)/Lambda{N_wav*M+1+Nt+j}(3); 0 0 0; 0 0 0];
    end
    
    Gamma{N_wav*M+1+Nt+Nt}(:,:,1)=[0; 0; 1];
    Gamma{N_wav*M+1+Nt+Nt}(:,:,2)=[0; 1; 0];
    Gamma{N_wav*M+1+Nt+Nt}(:,:,3)=[1; 0; 0];
     
   %% 
    for j=1:Nt
        j/Nt
       [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,N_wav*M+1+Nt+j,N_wav*M+1+2*j,maxSchmidtrank, Eps); 
    end
 
    %% SIMULATION

    [Gamma,Lambda,rhoS,siteS]...
        = FullEvolution(Gamma,Lambda,N_wav*M+1,UMPS,M, N_wav, Nt,maxSchmidtrank, Eps, showprogress);

    %% ANALYSIS

    exc_ev=zeros(3,Nt+1); 
    exc_ev(:,1)=abs(e_init).^2;
    [exc_ev]=GetExcitationProbabilities(rhoS,exc_ev);
    
    CorrelationTime=30;
    Ncorr=round(CorrelationTime/dt);

    a=diag(sqrt(1:dimB-1),1);
    autoc2=AutoCorrelation2(Gamma,Lambda,a',a,a',a,siteS-2*Ncorr-20,siteS-20);
    autoc1=AutoCorrelation1(Gamma,Lambda,a',a,length(Gamma)-2*Ncorr,length(Gamma)-2*Ncorr-2*Ncorr);
   
    autoc1=autoc1(1:2:end);
    
    autoc1=autoc1-real(autoc1(end));
    
    if showfigures==1
        figure
        plot((0:dt:TotalTime),exc_ev, '-')
        figure
        plot((0:dt:CorrelationTime),2*real(autoc2(1:2:end)), '-')
        figure
        plot((0:dt:CorrelationTime),real(autoc1), '-')
    end

     
     SpectrumLim=5;
     SpectrumRange=(-SpectrumLim:0.01:SpectrumLim);
     SpectrumV=zeros(1,length(SpectrumRange));
     
     for iter=1:length(SpectrumV)
         for iterM=1:length(autoc1)
             SpectrumV(iter)=SpectrumV(iter)+2*real(autoc1(iterM)*exp(1i*SpectrumRange(iter)*(iterM-1)*dt))*(TotalTime/gammap)^2/(64*pi);
         end
     end
     
     g2Function{itergammaP}=2*real(autoc2(1:2:end));
     SpectrumResult{itergammaP}=SpectrumV;

     
    if showfigures==1
        figure
        plot(SpectrumRange, SpectrumV);
        axis([-SpectrumLim SpectrumLim min(SpectrumV) max(SpectrumV)])
    end
    
end
if saveresults==1
    save([datestr(datetime('now')) '.mat'],'SpectrumResult','g2Function', 'gammapVec')
end