%% SETTINGS
clearvars
addpath('QSSE_functions/analysis/')
addpath('QSSE_functions/evolution/')
addpath('QSSE_functions/others/')
addpath('QSSE_functions/unitaries/')

OmegaVec=[5e-2];
SpectrumResult=cell(1,length(OmegaVec));
g2Function=cell(1,length(OmegaVec));

showprogress=0;
showfigures=0;
saveresults=1;

for iterOmega=1:length(OmegaVec)
    
OmegaL=OmegaVec(iterOmega);

TotalTime=120;        % Total simulation time
dt=0.05;          % Timestep

tau = 1;

N_wav=1;                                % Number of waveguides
N_sys=1;         % Number of separated components in the system
% x1=rand*exp(1i*rand*2*pi);
% x2=rand*exp(1i*rand*2*pi);
% e_init=[x1 x2 1]/sqrt(1+abs(x1)^2+abs(x2)^2);                         % Initial excitation of each state
e_init=[1 0 0];

dimB=3;                                 % Dimension for bosonic modes = 3 (no photon, 1 photon, 2 photons) for each waveguide
dimS=3;                                 % Dimension for each atom

Nt=round(TotalTime/dt);    % Total number of timesteps

maxSchmidtrank=120; % Maximum Schmidt rank for the SVD 
Eps = 10^(-12);      % Singular values truncation   

% Decay rates in the waveguide
gamma=1.;
gammaL=1*gamma; 
gammaR=1*gamma;

OmegaR=OmegaL;
myphi=0.1*2*pi;
    
DeltaL=-0;
DeltaR=0;
    

    m=cell(1,1);
    Theta=cell(1,1);
    UMPS=cell(1,1);
    
    % Number of steps between the interacting components
    m{1}=[-ceil(tau/dt) +1];

    % Losses associated
    Theta{1}=[0 0];
    
    %% UNITARY


    % Definition of Unitary

    UMPS{1}=V_system([gammaL gammaR], [OmegaL OmegaR], myphi, [DeltaL DeltaR], dimB,dt);

    %% INITIAL STATE

    m_max=1;
    N_loss=0;
    for j=1:length(m)
        m_max=max(m_max,max(max(abs(m{j}))));
        N_loss=N_loss+length(Theta{j}(Theta{j}~=0));
    end
    m_max=(N_wav+N_loss)*m_max;

    psiB=zeros(1,dimB);psiB(1)=1;           % Initial state of bosonic modes: vacuum
    Gammavac=reshape(psiB,[1,1,dimB]);
    Gammain=cell(1,N_wav);
    for j=1:N_wav+N_loss
        Gammain{j}=Gammavac;
    end

    GammaS=cell(1,N_sys);
    % Initial system state. [0,1]:excited
    for j=1:N_sys
        GammaS{j}=reshape(e_init,[1,1,dimS]);
    end


    Gamma0=cell(1,m_max);   % Initial Gamma
    Lambda0=cell(1,m_max);  % Initial Lambda
    for n=1:m_max
        Gamma0{n}=Gammavac;
        Lambda0{n}=1;
    end

    for j=1:N_sys
        Gamma0{m_max+j}=GammaS{j};
        Lambda0{m_max+j}=1;
    end

    %% SIMULATION

    [Gamma,Lambda,rhoS,rhoField, Error_Tracker, LambdaLoop]...
        = FullEvolution(Gamma0,Lambda0,Gammain,UMPS,m, N_wav,Theta,N_loss, Nt,maxSchmidtrank, Eps, showprogress);

    %% ANALYSIS
    exc_ev=zeros(3,Nt+1); 
    exc_ev(:,1)=abs(e_init).^2;
    [exc_ev]=GetExcitationProbabilities(rhoS,exc_ev);
    
    CorrelationTime=30;
    Ncorr=round(CorrelationTime/dt);

    a=diag(sqrt(1:dimB-1),1);
    autoc2=AutoCorrelation2(Gamma,Lambda,a',a,a',a,length(Gamma)-Ncorr-20,length(Gamma)-20);
    autoc1=AutoCorrelation1(Gamma,Lambda,a',a,length(Gamma)-20,length(Gamma)-Ncorr-20);
    
    autoc1=autoc1-(real(autoc1(end)));
    
    if showfigures==1
        figure
        plot((0:dt:TotalTime),exc_ev, '-')
        figure
        plot((0:dt:CorrelationTime),real(autoc2), '-')
        figure
        plot((0:dt:CorrelationTime),real(autoc1), '-')
    end

     
     SpectrumLim=4;
     SpectrumRange=(-SpectrumLim:pi/100:SpectrumLim);
     SpectrumV=zeros(1,length(SpectrumRange));
     
     for iter=1:length(SpectrumV)
         for iterM=1:length(autoc1)
             SpectrumV(iter)=SpectrumV(iter)+2*real(autoc1(iterM)*exp(1i*SpectrumRange(iter)*(iterM-1)*dt))/(OmegaL^4*2*pi);
         end
     end
     
     g2Function{iterOmega}=autoc2;
     SpectrumResult{iterOmega}=SpectrumV;
     
    if showfigures==1
        figure
        plot(SpectrumRange, SpectrumV);
        axis([-SpectrumLim SpectrumLim min(SpectrumV) max(SpectrumV)])
    end
     
end

if saveresults==1
    save([datestr(datetime('now')) '.mat'],'SpectrumResult','g2Function')
end

%%
%  figure
%  plot((0:dt:Totaltime),exc_ev) 
%  legend('g population','r population','l population')