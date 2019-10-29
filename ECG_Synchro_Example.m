%%--------------------   ECG DATA ANALYSIS  ----------------------
%-----------------------------------------------------------------
% 08/02/2019 : Verdière. K. : Creation
%-----------------------------------------------------------------

% Script for synchronization analysis of cardiac data (R peaks).
% (See Verdiere, K. J., Albert, M., Dehais, F., and Roy, R. N. (2019). Physiological
% synchrony revealed by delayed coincidence count: Application to a cooperative complex
% environment. IEEE Transactions on Human-Machine Systems.)
% User feedback welcome, email: kevin@verdiere.fr
% 
% The permutation method is from 
% Albert, M. (2015). Tests of independence by bootstrap and permutation: an asymptotic
% and non-asymptotic study. Application to neurosciences. PhD thesis, Université Nice Sophia
% Antipolis.
% 
% Authors: Kevin Verdière ISAE-SUPAERO Toulouse France
%
% Description:
%   This script computes an example synchronization analysis for three couples
%   of Pilot Monitoring (PM) and Pilot Flying (PF). 
% 
% Usage: Launch the entire script, make sure RR_example.mat is in the
% script folder. 

clear
load('Rpeak_example.mat')
% The R variable is a cell contains the detected R peaks.  
% R{Subject, Scenario, Teammates}
% ----- Subject: representes the couple, here we have three couples of two
% participants
% ----- Scenario: represents the experimental conditions. Scenario are
% evaluated indeptendently (here we have 8 conditions).
% ----- Teammates: A team is composed of 2 teammate: One PF, one PM. (must be 2)
%% Synchronization
% COICIDENCE CALCULATION 

N_SujMax    = size(R,1); % Number of subjects
N_Scenario  = size(R,2); % Number of scenario
Delta       = 20; % Time in millisecond in which R peak are considerend coincident.
N_Perm      = 100000; % Number of permutation


Coincidence=nan(N_SujMax,N_SujMax,N_Scenario); % Init Variable
t1=datetime;
for Scenario=1:N_Scenario
    
    IncremSuj=0;
    for CoupPF=1:N_SujMax
        for CoupPM=1:N_SujMax
            clc; IncremSuj=IncremSuj+1;
            disp(['Coincidence Calculation...'])
            disp(['Scenario ' num2str(Scenario) ' / ' num2str(N_Scenario)])
            disp(['Couple PF ' num2str(CoupPF) ' <-> PM ' num2str(CoupPM)])
            
            xPF=round(R{CoupPF,Scenario,1}*1000); % convert to ms
            xPM=round(R{CoupPM,Scenario,2}*1000); % convert to ms
            
            if ~isempty(xPM) && ~isempty(xPF)
                PF=xPF';
                PM=xPM;
                
                PF=repmat(PF,1,size(PM,2));
                PM=repmat(PM,size(PF,1),1);
                
                Coincidence(CoupPF,CoupPM,Scenario)=length(find(abs(PF-PM)<Delta));
            end
            
        end
    end
end
t2=datetime;
%% PERMUTATION
datapvalue=[];
tb=[];

CoupleUsed=1:N_SujMax; % Select participants you want to include. Here all
N_CoupleUsed=length(CoupleUsed); 

t3=datetime;

%----- Find if there is NaN in data, and supress participants if so
for i=1:N_Scenario
    Scenar=Coincidence(CoupleUsed,CoupleUsed,i);
    
    PFdead=find(all(isnan(Scenar),2));
    PMdead=find(all(isnan(Scenar)));
    
    MaskSuj=1:N_CoupleUsed;
    MaskSuj(sort([PFdead ; PMdead']))=[];
    
    Data{i}=Scenar(MaskSuj,MaskSuj);
    
end

%----- Permutation
disp(['Permutation : ' num2str(N_Perm)])
Result=nan(N_Perm,N_Scenario);
Diago=[];

for S=1:N_Scenario
    disp(['Scenario ' num2str(S)])
    NSuj=size(Data{S},1);
    EyeMat=logical(eye(NSuj));
    Diago(S)=trace(Data{S});
    
    for i=1:N_Perm
        Mask=EyeMat(randperm(NSuj),:);
        Result(i,S)=sum(Data{S}(Mask));
    end
    
end
Diago=repmat(Diago,N_Perm,1);

% pplus and pmoins contains raw p-value. 
% if pplus  is low (<.05 for example) there is a SIGNIFICANT SYNCHRONIZATION
% if pmoins is low (<.05 for example) there is a SIGNIFICANT DESYNCHRONIZATION
pplus   = (sum(Result>=Diago)+1) / (N_Perm+1);
pmoins  = (sum(Result<=Diago)+1) / (N_Perm+1);

datapvalue(:,:)=[pmoins ; pplus]';

[Val{1},Ordre{1}]=sort(pplus);
[Val{2},Ordre{2}]=sort(pmoins);

NbSignif = 0;
for SynchDesynch = 1:2
    tmpVal = Val{SynchDesynch};
    
    tb(:,1:2)=[tmpVal; (1:length(tmpVal))*0.05/length(tmpVal)]';
    tb(:,3)=tb(:,1)<=tb(:,2);
    tb(:,4)=Ordre{SynchDesynch};
    Results{SynchDesynch} = array2table(tb,...
        'VariableNames',{'pValue','BenjaminiThreshold','UnderThreshold','ScenarioNumber'});
    
    NbSignif = NbSignif + length(find(tb(:,3)));
end
t4=datetime;

Results{1} % Synchro results
Results{2} % Desynchro results

disp(['There is ' num2str(NbSignif) ' Scenario significantly different' ...
    '(i.e. more or less synchronized) for real teams compared to '...
    'permutated teams'])


%% Graph representing each scenario
figure
b=bar([Diago(1,:) ; mean(Result)]');
grid on; hold on

set(gcf,'position',[12  194  979  722])
errorbar([1:N_Scenario]+.15, mean(Result) , std(Result) ,'xk')

legend({'C^{obs} (Real Teams)','C^{b} (Permutated Teams)'},'location','southeast')
xlabel('Scenarii'); ylabel('Total Coincidence Count (trace)')


set(gca,'fontsize',18)
b(1).FaceColor=[1 1 1]*.2;
b(2).FaceColor=[1 1 1]*.8;



