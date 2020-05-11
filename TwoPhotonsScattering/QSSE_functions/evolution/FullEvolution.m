function [Gamma,Lambda,rhoS,siteS]...
            =FullEvolution(Gamma,Lambda,siteS,UMPS,M, N_wav, Nt,maxSchmidtrank, Eps, showprogress)
% [Gamma,Lambda,rhoS, rhoField, ErrorTracker, LambdaLoop]=FullEvolution(Gamma,Lambda,Gammain,UMPS,m, N_wav,Theta, N_loss,Nt,maxSchmidtrank, Eps, showprogress)
% 
% Evolves the Gamma and Lambda tensors according to the unitary UMPS for Nt
% timesteps. 
%
% Gammamin is a cell which contains the initial states of every
% new photonic bin. m is the cell array containing the matrices for the time-delay of the interacting time-bins.
% Theta is the loss cell array, and N_loss the
% number of loss channels. 
%
% If showprogress is 1, displays a waitbar showing
% the progress of the algorithm, and prints the computation time.
%
% maxSchmidtrank is the maximal authorised number of singular values in the
% Schmidt decomposition. Eps is the minimal value for each singular
% value considered non-zero

warning('on','mytest:maxrank')

rhoS=cell(1,Nt);

%% Time evolution

if showprogress
    prog_bar = waitbar(0);
    titleHandle = get(findobj(prog_bar,'Type','axes'),'Title');
    set(titleHandle,'FontSize',20);
    waitbar(0,prog_bar,sprintf('0.00%',0));
    tic
end

for k=1:Nt 

           
    [Gamma,Lambda]=...
        Shift_site1_to_siteN(Gamma,Lambda,siteS-N_wav*M,siteS-1,maxSchmidtrank,Eps);
    
    [Gamma,Lambda]=...
        Shift_site1_to_siteN(Gamma,Lambda,siteS,siteS+2,maxSchmidtrank,Eps);

    % Apply the local Unitary evolution
   [Gamma,Lambda]=LocalUnitaryEvolution(Gamma,Lambda,siteS+2,UMPS,maxSchmidtrank,Eps);

    % State analysis
    rhoS{1,k}=single_site_reduced_state_Gamma_efficient(Gamma,Lambda,siteS+2);

    [Gamma,Lambda]=...
        Shift_site1_to_siteN(Gamma,Lambda,siteS-1,siteS-N_wav*M,maxSchmidtrank,Eps);
            
     siteS=siteS+2;% Update system index
         
    % Display progression
    if showprogress
        waitbar(k/Nt,prog_bar,sprintf('%3.2f%%',100*k/Nt));
    end
end

if showprogress
    toc
    close(prog_bar)
end