function [Gamma,Lambda,rhoS, rhoField, ErrorTracker, LambdaLoop]...
            =FullEvolution(Gamma,Lambda,Gammain,UMPS,m, N_wav,Theta, N_loss,Nt,maxSchmidtrank, Eps, showprogress)
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
N_sys=length(m);

%% Checking the correctness of matrix m

check_order=zeros(1,N_wav);
for j=1:N_sys
    if find(m{j}>1)
       error('The elements of the delay matrix M should be at most +1') 
    elseif size(m{j},1)~= N_wav
       error('The delay matrix M should contain N_wav rows (eventually filled with 0), excluding the losses waveguides') 
    elseif size(Theta{j})~=size(m{j})
        error('The loss matrix THETA and the delay matrix M should be of the same size')
    elseif length(size(UMPS{j}))~= 2*(length(m{j}(m{j}~=0))+1)
        error('The number of indices of UMPS is not consistent with the number of non-zero entries in the delay matrix M')
    end
    for i=1:N_wav
        if ~issorted(m{j}(i,:)) || length(unique(find(m{j}(i,:))))<length(find(m{j}(i,:)))
            error('Each row of the delay matrix M should be sorted in strictly ascending order')
        elseif find(m{j}(i,:)==1)
            if check_order(i) ~=0
                error('There should be at most one "+1" per waveguide in the delay matrix M')
            else check_order(i)=j;
            end
        end
    end
    for j2=j+1:N_sys
        for i=1:N_wav
           if numel(intersect(m{j}(i,:),m{j2}(i,:)))~=0
               error('The rows of the delay matrix m{j} should not contain any common element with those of matrix m{j2} for j~=j2')
           end
        end
    end
end
if find(check_order==0)
    error('There should be at least one "+1" per waveguide in the delay matrix M')
end

%% Reordering the system components to match the waveguide order

siteS=length(Gamma)-(N_sys-1); % Index for the system

order_index=check_order;

for j=1:N_wav
    if order_index(j)>j
           [Gamma,Lambda]=Shift_site1_to_siteN(Gamma,Lambda,siteS+order_index(j)-1,siteS+(j-1),maxSchmidtrank,Eps);
           order_index(order_index==order_index(j))=j;
           temp=order_index(j+1:end);
           temp(temp<order_index(j))=temp(temp<order_index(j))+1;
           order_index(j+1:end)=temp;
    end
end

vec_sys = (1:N_sys);
vec_sys(check_order)=[];
order_index=[check_order vec_sys];

%% Defining the sites interacting with each system component

siteTimeBinsL=cell(1,N_sys);
siteTimeBinsR=cell(1,N_sys);
temp_mat=cell(1,N_sys);
minSite=siteS-1;

for j=1:N_sys
    temp_mat{j}=m{order_index(j)}*(N_wav+N_loss)+repmat((0:N_wav-1)',1,size(m{order_index(j)},2));
    siteTimeBinsL{j}=siteS+temp_mat{j}(temp_mat{j}<0);
    siteTimeBinsR{j}=siteS+temp_mat{j}(temp_mat{j}>(N_wav+N_loss)-1)-(N_wav+N_loss)+N_sys;
    minSite=min(minSite,min(siteTimeBinsL{j}));
end

rhoS=cell(N_sys,Nt);
rhoField=cell(N_sys,length(siteTimeBinsL{j})+length(siteTimeBinsR{j}),Nt);

%% Defining the bins which undergo loss

site_loss=cell(1,N_sys);
theta_loss=cell(1,N_sys);
time_loss=cell(1,N_sys);

for j=1:N_sys
    site_loss{j}=temp_mat{j}(Theta{order_index(j)}~=0);
    site_loss{j}(site_loss{j}>0)=site_loss{j}(site_loss{j}>0)-(N_wav+N_loss);
    site_loss{j}=siteS+site_loss{j};
    
    theta_loss{j}=Theta{order_index(j)}(Theta{order_index(j)}~=0);
    
    time_loss{j}=m{order_index(j)};
    time_loss{j}(time_loss{j}>0)=time_loss{j}(time_loss{j}>0)-1;
    time_loss{j}=siteS+time_loss{j}*(N_wav+N_loss);
end

U_loss=cell(N_sys,1);
for j=1:N_sys    
    for l=1:length(site_loss{j})
        U_loss{j,l}=U_WaveguideLosses(size(Gammain{1},3),theta_loss{j}(l));
    end
end

%% Time evolution

if showprogress
    prog_bar = waitbar(0);
    titleHandle = get(findobj(prog_bar,'Type','axes'),'Title');
    set(titleHandle,'FontSize',20);
    waitbar(0,prog_bar,sprintf('0.00%',0));
    tic
end

LambdaLoop=cell(1,Nt);
ErrorTracker=ones(Nt,2);
for k=1:Nt 
    
     for j=1:N_wav+N_loss %New empty bins
        Gamma{siteS+j+(N_sys-1)}=Gammain{j};
        Lambda{siteS+j+(N_sys-1)}=1;
     end

     for j=1:N_sys %Each subsystem evolves
         
            for l=length(siteTimeBinsL{j}):-1:1 % SWAP the feedback bins close to the system
                [Gamma,Lambda,ErrorT]=...
                    Shift_site1_to_siteN(Gamma,Lambda,siteTimeBinsL{j}(l),siteS-length(siteTimeBinsL{j})+(l-1),maxSchmidtrank,Eps);
                ErrorTracker(k,:)=[ErrorTracker(k,1)*ErrorT(1),min(ErrorTracker(k,2),ErrorT(2))];
            end
            for l=1:length(siteTimeBinsR{j}) % SWAP the feedback bins close to the system
                [Gamma,Lambda,ErrorT]=...
                    Shift_site1_to_siteN(Gamma,Lambda,siteTimeBinsR{j}(l),siteS+(l-1),maxSchmidtrank,Eps);
                ErrorTracker(k,:)=[ErrorTracker(k,1)*ErrorT(1),min(ErrorTracker(k,2),ErrorT(2))];
            end

            % Apply the local Unitary evolution
           [Gamma,Lambda,ErrorT]=LocalUnitaryEvolution(Gamma,Lambda,siteS+length(siteTimeBinsR{j}),UMPS{order_index(j)},maxSchmidtrank,Eps);
                ErrorTracker(k,:)=[ErrorTracker(k,1)*ErrorT(1),min(ErrorTracker(k,2),ErrorT(2))];
            
            % State analysis
            rhoS{order_index(j),k}=single_site_reduced_state_Gamma_efficient(Gamma,Lambda,siteS+length(siteTimeBinsR{j}));

            for l=1:length(siteTimeBinsL{j}) % SWAP back the bins to their original positions
                [Gamma,Lambda,ErrorT]=...
                    Shift_site1_to_siteN(Gamma,Lambda,siteS-length(siteTimeBinsL{j})+(l-1),siteTimeBinsL{j}(l),maxSchmidtrank,Eps);
                ErrorTracker(k,:)=[ErrorTracker(k,1)*ErrorT(1),min(ErrorTracker(k,2),ErrorT(2))];
            end
            for l=1:length(siteTimeBinsR{j}) % SWAP back the bins to their original positions
                [Gamma,Lambda,ErrorT]=...
                    Shift_site1_to_siteN(Gamma,Lambda,siteS+(l-1),siteS+(l-1)-(j-1),maxSchmidtrank,Eps);
                ErrorTracker(k,:)=[ErrorTracker(k,1)*ErrorT(1),min(ErrorTracker(k,2),ErrorT(2))];
            end
            
             siteS=siteS+length(siteTimeBinsR{j})+1;% Update system index
             
     end
     
    siteS=siteS-N_sys;

    for l=1:N_loss
        [Gamma,Lambda,ErrorT]=...
            Shift_site1_to_siteN(Gamma,Lambda,siteS+N_sys+(l-1),siteS+(l-1),maxSchmidtrank,Eps);
        ErrorTracker(k,:)=[ErrorTracker(k,1)*ErrorT(1),min(ErrorTracker(k,2),ErrorT(2))];
    end
    
    siteS=siteS+N_loss;
    
    for j=1:N_sys
        for l=1:length(site_loss{j})
            [Gamma,Lambda,ErrorT]=...
                Shift_site1_to_siteN(Gamma,Lambda,time_loss{j}(l)+N_wav,site_loss{j}(l)+1,maxSchmidtrank,Eps);
            ErrorTracker(k,:)=[ErrorTracker(k,1)*ErrorT(1),min(ErrorTracker(k,2),ErrorT(2))];
            
            [Gamma,Lambda,ErrorT]=...
                LocalUnitaryEvolution(Gamma,Lambda,site_loss{j}(l)+1,U_loss{j,l},maxSchmidtrank,Eps);
            ErrorTracker(k,:)=[ErrorTracker(k,1)*ErrorT(1),min(ErrorTracker(k,2),ErrorT(2))];

            [Gamma,Lambda,ErrorT]=...
                Shift_site1_to_siteN(Gamma,Lambda,site_loss{j}(l)+1,time_loss{j}(l)+N_wav+N_loss-1,maxSchmidtrank,Eps);
            ErrorTracker(k,:)=[ErrorTracker(k,1)*ErrorT(1),min(ErrorTracker(k,2),ErrorT(2))];
        end  
    end
    
    Lambda_temp=Lambda(minSite:siteS-1);
    LambdaLoop{k}=EntropyProfile(Lambda_temp);
    
    for j=1:N_sys
        for l=1:length(siteTimeBinsL{j})
            rhoField{order_index(j),l,k}=single_site_reduced_state_Gamma_efficient(Gamma,Lambda,siteTimeBinsL{j}(l));
        end
        for l=1:length(siteTimeBinsR{j})
            rhoField{order_index(j),length(siteTimeBinsL{j})+l,k}=single_site_reduced_state_Gamma_efficient(Gamma,Lambda,siteTimeBinsR{j}(l)-N_sys);
        end
    end
    
    for j=1:N_sys
         siteTimeBinsR{j}=siteTimeBinsR{j}+N_wav+N_loss;
         siteTimeBinsL{j}=siteTimeBinsL{j}+N_wav+N_loss;
         site_loss{j}=site_loss{j}+N_wav+N_loss;
         time_loss{j}=time_loss{j}+N_wav+N_loss;
    end
    
    minSite=minSite+N_wav+N_loss;
         
    % Display progression
    if showprogress
        waitbar(k/Nt,prog_bar,sprintf('%3.2f%%',100*k/Nt));
    end
end

if showprogress
    toc
    close(prog_bar)
end