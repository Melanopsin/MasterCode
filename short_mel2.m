% Melanopsin(totrun,dataset) will run the Gillespie algorithm code
% for the number of runs specified (totrun) and using the initial values 
% and deterministic rate constants for the specified dataset (dataset). Later once
% the optimization code is finished the rate constants used will be from
% the rates written to file by that routine. The realizations are saved in files called
% run<number>.mat.

% % function short_mel(totrun,dataset,tmax,flashint)
% tmax: final time for each run
clear all
whatha=0;
dataset=5
%% I move the location of the code which save the parameter mat file
%% so that the file does not save additional variables.
switch dataset    
    
    case 1
        Set_IC_ephys
        load('ephys.mat')     % electrophysiology
        Edata = dlmread('EphysGraph.csv');
        exp_data = [Edata(1100:end,2)-1.1 Edata(1100:end,5)];

    case 2
        Set_ICMulti
        load('multidata.mat') % multi flash
    case 3
        Set_IC_beta1
        load('beta1.mat')     % overexpressed beta arrestin1
    case 4
        Set_IC_beta2
        load('beta2.mat')     % overexpressed beta arrestin2 
    case 5 
        Set_IC_clhkah2
        load('clhkah.mat')   % calcium imaging is the default    
    otherwise 
        Set_IC_hwk
        load('hwk.mat')   % calcium imaging is the default
        
end

totrun=1000
tmax=15
flashint=100

% % clf
%% to set a random seed
seed_name = sum(100*clock);
rand('seed',seed_name)
%% to set a random seed

%% store time, molecule numbers in every 'time_step' sec for all runs
time_step=0.0625; % maximum time step, before we set as 0.25
counter =1; % counter counts the number of time iterations.
maxcounter=10000000;

%% You record t, M, X in every time interval equal to time_step
Xstore = zeros(floor(tmax/time_step)+1,12,totrun); 
Mstore = zeros(floor(tmax/time_step)+1,49,totrun); 
ttstore = zeros(floor(tmax/time_step)+1,1,totrun);

no_rxns = 88;                   % (<--- this is hardcode) number of reactions (total)
reaction_count = zeros(no_rxns,totrun); % it counts the number of each reactions' occuring.
%% store time, molecule numbers in every 'time_step' sec for all runs
slowmedown=1.00;
tic;

%% start of 'totrun' simulations
for runnum = 1:totrun

max_chan = 0;       % What are these?
maxchan_time = 0;   % What are these?

t = 0;
tflash = flashint;                   % Time for second flash
tjump = tflash;                 % Time between flashes

% MODEL OF MELANOPSIN ACTIVATION

% Mn* + G.GDP <-- kG1 y(n)/kG2 --> Mn*.G.GDP
% Mn*.G.GDP -- kG3 --> Mn*.G + GDP
% Mn*.G + GTP -- kG4 --> Mn*.G.GTP
% Mn*.G.GTP -- kG5 --> Mn* + Ga.GTP + Gbg
% PLC + Ga.GTP -- kP --> PLC*.Ga.GTP
% PLC*.Ga.GTP -- kI1 --> PLC.Ga.GDP
% PLC.Ga.GDP -- kI2 --> PLC + Ga.GDP
% Ga.GDP + Gbr -- kI3 --> G.GDP
% PIP2 + PLC*.Ga.GTP -- kS --> SecM + PLC*.Ga.GTP
% SecM -- delta --> 0
% SecM + Channel+ <-- kO/kC --> SecM.Channel-
% Mn* + K <-- kK1/kK2 --> Mn*.K
% Mn*.K + ATP -- kK3 --> Mn+1*.K + ADP
% Mn* + ArrB1 -- kB1 w(n) --> Mn.ArrB1
% Mn* + ArrB2 -- kB2 w(n) --> Mn.ArrB2
% Mn.ArrB1 -- kUB1 w(n) --> Mp + ArrB1
% Mn.ArrB2 -- kUB2 w(n) --> Mp + ArrB2
% Mp -- kDe w(6) --> M0


%% initialize the parameters and initial conditions in each run 
switch dataset    
    
    case 1
        load('ephys.mat')     % electrophysiology
    case 2
        load('multidata.mat') % multi flash
    case 3
        load('beta1.mat')     % overexpressed beta arrestin1
    case 4
        load('beta2.mat')     % overexpressed beta arrestin2  
    case 5
        load('clhkah.mat')
    otherwise 
        load('hwk.mat')   % calcium imaging is the default
end


%% SPECIES: X = [
%% X(1)           G.GDP
%% X(2)           Ga.GTP
%% X(3)           Gbg
%% X(4)           PLC
%% X(5)           PLC*.Ga.GTP
%% X(6)           SecM      
%% X(7)           Channel-
%% X(8)           SecM.Channel+
%% X(9)           PLC.Ga.GDP
%% X(10)          Ga.GDP    ];
%% X(11)          Beta1;
%% X(12)          Beta2;


%% MELANOPSIN COMPLEXES: M = [
%% M(1)                        M0*
%% M(2)                        M0*.G.GDP
%% M(3)                        M0*.G
%% M(4)                        M0*.G.GTP
%% M(5)                        M0*.K
%% M(6)                        M1*
%% M(7)                        M1*.G.GDP
%% M(8)                        M1*.G
%% M(9)                        M1*.G.GTP
%% M(10)                       M1*.K
%% M(11)                       M1*.ArrB1
%% M(12)                       M1*.ArrB2
%% M(13)                       M2*
%% M(14)                       M2*.G.GDP
%% M(15)                       M2*.G
%% M(16)                       M2*.G.GTP
%% M(17)                       M2*.K
%% M(18)                       M2*.ArrB1
%% M(19)                       M2*.ArrB2
%% M(20)                       M3*
%% M(21)                       M3*.G.GDP
%% M(22)                       M3*.G
%% M(23)                       M3*.G.GTP
%% M(24)                       M3*.K
%% M(25)                       M3*.ArrB1
%% M(26)                       M3*.ArrB2
%% M(27)                       M4*
%% M(28)                       M4*.G.GDP
%% M(29)                       M4*.G
%% M(30)                       M4*.G.GTP
%% M(31)                       M4*.K
%% M(32)                       M4*.ArrB1
%% M(33)                       M4*.ArrB2
%% M(34)                       M5*
%% M(35)                       M5*.G.GDP
%% M(36)                       M5*.G
%% M(37)                       M5*.G.GTP
%% M(38)                       M5*.K
%% M(39)                       M5*.ArrB1
%% M(40)                       M5*.ArrB2
%% M(41)                       M6*
%% M(42)                       M6*.G.GDP
%% M(43)                       M6*.G
%% M(44)                       M6*.G.GTP
%% M(45)                       M6*.K
%% M(46)                       M6*.ArrB1
%% M(47)                       M6*.ArrB2 ];
%% M(48)                       M0
%% M(49)                       MP



K = [ kG1,      kG2,        kG3,        kG4*GTP,    kG5, ...    
      kP,       kI1,        kS*PIP2,    kO,         kC, ...
      kk1*Ki,   kk2,        kk3*ATP,    kB1,        kB2, ...
      kmax,     KM,         kI2,        kI3, ... 
      kUB1      kUB2        kDe];
  
%%  K(1)=kG1,       K(2)=kG2,     K(3)=kG3,      K(4)=kG4*GTP,      K(5)=kG5, 
%%  K(6)=kP,        K(7)=kI1,     K(8)=kS*PIP2,  K(9)=kO,           K(10)=kC,
%%  K(11)=kk1*Ki,   K(12)=kk2,    K(13)=kk3*ATP, K(14)=kB1,   K(15)=kB2,  
%%  K(16)=kmax,     K(17)=KM,     K(18)=kI2,     K(19)=kI3         
%%  K(20) = kUB1    K(21) = kUB2  K(22) = kDe


%% set initial values
h_tot=0; % initialize h_tot
h = zeros(no_rxns,1);  % initialize the hazard vector
tstore(1,1,runnum) = t;
Xstore(1,:,runnum) = X;
Mstore(1,:,runnum) = M;
ttstore(1,1,runnum) = 0;
prev_t_index = 1; % it stores the previous time intex
%% set initial values
runnum

%% Build function to account for increase in arrestin binding affinity
%% with more phosphates bound to melanopsin carboxyl tail.

% W = @(n) 1-exp(-n*10);
W = @(n) 1-exp(-n*100);

%% Build function to account for decrease in G-protein activation
%% with more phosphates bound to melanopsin carboxyl tail.

% Y = @(n) exp(-n);
Y = @(n) exp(-n/1000);
Z = @(n) exp(-2*n); % for decreasing phosphorylation rate kK1

%% Begin the algorithm

for counter=1:maxcounter
    % check if the final time has been reached or exceeded
    if t>=tmax
        break;
    end
    t;
     if M(1)<10
         slowmedown=1.00;
     end
    
     if t>= tflash
        % signals the next flash
       tflash = tflash + tjump;
       % the signal hits all of the melanopsin, potentially activating 10%
       % at any given flash. For a first approximation, assume that the 10%
       % does distinguish between activated and inactive...Pool the M0 and
       % M0*
       numactivated = floor((M(48))*0.1); % For now just activate 10% every time floor it to keep from going negative   
       M(48) = M(48)-numactivated; % do it this way to keep it an integer and not risk losing 
       if M(48) < 0
           M(48) = 0
       else
       % M0s from floors and ceilings.
       M(1) = M(1) + numactivated; % same thing
       end
    end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % use cell #'s to construct hazards for the current time step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% G-Protein Activation
    
    % M0* + G.GDP -- kG1*y(n) --> M0*.G.GDP
    h(1)  = Y(0)*K(1)*M(1)*X(1);    
    % And now for the other phosphorylated species
    h(6)  = Y(1)*K(1)*M(6)*X(1);
    h(11) = Y(2)*K(1)*M(13)*X(1);
    h(16) = Y(3)*K(1)*M(20)*X(1);
    h(21) = Y(4)*K(1)*M(27)*X(1);
    h(26) = Y(5)*K(1)*M(34)*X(1);
    h(31) = Y(6)*K(1)*M(41)*X(1);

    % M0*.G.GDP -- kG2 --> M0* + G.GDP
    h(2)  = K(2)*M(2);         
    % And now for the other phosphorylated species
    h(7)  = K(2)*M(7);
    h(12) = K(2)*M(14);
    h(17) = K(2)*M(21);
    h(22) = 0*K(2)*M(28);
    h(27) = 0*K(2)*M(35);
    h(32) = 0*K(2)*M(42);

    % M0*.G.GDP -- kG3 --> M0*.G + GDP
    h(3)  = K(3)*M(2);  
    % And now for the other phosphorylated species
    h(8)  = K(3)*M(7);
    h(13) = K(3)*M(14);
    h(18) = K(3)*M(21);
    h(23) = 0*K(3)*M(28);
    h(28) = 0*K(3)*M(35);
    h(33) = 0*K(3)*M(42);

    % M0*.G + GTP -- kG4 --> M0*.G.GTP
    h(4)  = K(4)*M(3);        
    % And now for the other phosphorylated species
    h(9)  = K(4)*M(8);
    h(14) = K(4)*M(15);
    h(19) = K(4)*M(22);
    h(24) = 0*K(4)*M(29);
    h(29) = 0*K(4)*M(36);
    h(34) = 0*K(4)*M(43);

    % M0*.G.GTP -- kG5 --> M0* + Ga.GTP + Gbg
    h(5)  = K(5)*M(4);         
    % And now for the other phosphorylated species
    h(10) = K(5)*M(9);
    h(15) = K(5)*M(16);
    h(20) = K(5)*M(23);
    h(25) = 0*K(5)*M(30);
    h(30) = 0*K(5)*M(37);
    h(35) = 0*K(5)*M(44);
    
%%
% PLC and G-protein activation/inactivation
    
    h(36) = K(6)*X(4)*X(2);         % PLC + Ga.GTP -- kP --> PLC*.Ga.GTP
    h(37) = K(7)*X(5);              % PLC*.Ga.GTP -- kI1 --> PLC.Ga.GDP
    h(74) = K(18)*X(9);             % PLC.Ga.GDP -- kI2 --> PLC + Ga.GDP
    h(75) = K(19)*X(10)*X(3);       % Ga.GDP + Gbr -- kI3 --> G.GDP
                                    
%%
% Second Messenger Creation
    
    h(38) = K(8)*X(5);              % PIP2 + PLC*.Ga.GTP -- kS --> SecM + PLC*.Ga.GTP
    
%%
% Channel Opening/Closing
    
    h(39) = K(9)*X(6)*X(7);         % SecM + Channel- -- kO --> SecM.Channel+
    h(40) = K(10)*X(8);             % SecM.Channel+ -- kC --> SecM + Channel-
    
%%
% Kinase Phosphorylation
    
    h(41) = Z(0)*K(11)*M(1);          % M0* + K -- kk1 --> M0*K 
    h(42) = K(12)*M(5);               % M0*K -- kk2 --> M0* + K
    h(43) = K(13)*M(5);               % M0*K + ATP -- kk3 --> M1*K + ADP
    
    h(44) = Z(1)*K(11)*M(6);          % M1* + K -- kk1 --> M1*K 
    h(45) = K(12)*M(10);              % M1*K -- kk2 --> M1* + K
    h(46) = K(13)*M(10);              % M1*K + ATP -- kk3 --> M2*K + ADP
    
    h(47) = Z(2)*K(11)*M(13);         % M2* + K -- kk1 --> M2*K 
    h(48) = K(12)*M(17);              % M2*K -- kk2 --> M2* + K
    h(49) = K(13)*M(17);            % M2*K + ATP -- kk3 --> M3*K + ADP
    
    h(50) = Z(3)*K(11)*M(20);       % M3* + K -- kk1 --> M3*K 
    h(51) = K(12)*M(24);            % M3*K -- kk2 --> M3* + K
    h(52) = 0*K(13)*M(24);            % M3*K + ATP -- kk3 --> M4*K + ADP
    
    h(53) = 0*Z(4)*K(11)*M(27);       % M4* + K -- kk1 --> M4*K 
    h(54) = 0*K(12)*M(31);            % M4*K -- kk2 --> M4* + K
    h(55) = 0*K(13)*M(31);            % M4*K + ATP -- kk3 --> M5*K + ADP
    
    h(56) = 0*Z(5)*K(11)*M(34);       % M5* + K -- kk1 --> M5*K 
    h(57) = 0*K(12)*M(38);            % M5*K -- kk2 --> M5* + K
    h(58) = 0*K(13)*M(38);            % M5*K + ATP -- kk3 --> M6*K + ADP
            
    h(59) = 0*Z(6)*K(11)*M(41);       % M6* + K -- kk1 --> M6*K 
    h(60) = 0*K(12)*M(45);            % M6*K -- kk2 --> M6* + K
    
%%
% Arrestin Binding
    
    h(61) = W(1)*K(14)*M(6)*X(11)*slowmedown;        % M1* + ArrB1 -- kB1*w(n) --> M1.ArrB1
    h(62) = W(1)*K(15)*M(6)*X(12)*slowmedown;        % M1* + ArrB2 -- kB2*w(n) --> M1.ArrB2
    
    h(63) = W(2)*K(14)*M(13)*X(11)*slowmedown;       % M2* + ArrB1 -- kB1*w(n) --> M2.ArrB1
    h(64) = W(2)*K(15)*M(13)*X(12)*slowmedown;       % M2* + ArrB2 -- kB2*w(n) --> M2.ArrB2

    h(65) = W(3)*K(14)*M(20)*X(11)*slowmedown;       % M3* + ArrB1 -- kB1*w(n) --> M3.ArrB1
    h(66) = W(3)*K(15)*M(20)*X(12)*slowmedown;       % M3* + ArrB2 -- kB2*w(n) --> M3.ArrB2
    
    h(67) = 0*W(4)*K(14)*M(27)*X(11)*slowmedown;       % M4* + ArrB1 -- kB1*w(n) --> M4.ArrB1
    h(68) = 0*W(4)*K(15)*M(27)*X(12)*slowmedown;       % M4* + ArrB2 -- kB2*w(n) --> M4.ArrB2
    
    h(69) = 0*W(5)*K(14)*M(34)*X(11)*slowmedown;       % M5* + ArrB1 -- kB1*w(n) --> M5.ArrB1
    h(70) = 0*W(5)*K(15)*M(34)*X(12)*slowmedown;       % M5* + ArrB2 -- kB2*w(n) --> M5.ArrB2
    
    h(71) = 0*W(6)*K(14)*M(41)*X(11)*slowmedown;       % M6* + ArrB1 -- kB1*w(n) --> M6.ArrB1
    h(72) = 0*W(6)*K(15)*M(41)*X(12)*slowmedown;       % M6* + ArrB2 -- kB2*w(n) --> M6.ArrB2
    
%%
% SecM degradation
    h(73) = kmax*X(6)/(X(6)+KM);    % SecM -- delta --> 0, delta = kmax*SecM/(SecM+KM)
    
%%
% Arrestin unbinding and deactivation

    h(76) = W(1)*K(20)*M(11);       % M1.ArrB1 -- kUB1 --> MP + ArrB1
    h(77) = W(2)*K(20)*M(18);       % M2.ArrB1 -- kUB1 --> MP + ArrB1
    h(78) = W(3)*K(20)*M(25);       % M3.ArrB1 -- kUB1 --> MP + ArrB1
    h(79) = 0*W(4)*K(20)*M(32);       % M4.ArrB1 -- kUB1 --> MP + ArrB1
    h(80) = 0*W(5)*K(20)*M(39);       % M5.ArrB1 -- kUB1 --> MP + ArrB1
    h(81) = 0*W(6)*K(20)*M(46);       % M6.ArrB1 -- kUB1 --> MP + ArrB1
    
    h(82) = W(1)*K(21)*M(12);       % M1.ArrB2 -- kUB2 --> MP + ArrB2
    h(83) = W(2)*K(21)*M(19);       % M2.ArrB2 -- kUB2 --> MP + ArrB2
    h(84) = W(3)*K(21)*M(26);       % M3.ArrB2 -- kUB2 --> MP + ArrB2
    h(85) = 0*W(4)*K(21)*M(33);       % M4.ArrB2 -- kUB2 --> MP + ArrB2
    h(86) = 0*W(5)*K(21)*M(40);       % M5.ArrB2 -- kUB2 --> MP + ArrB2
    h(87) = 0*W(6)*K(21)*M(47);       % M6.ArrB2 -- kUB2 --> MP + ArrB2
    
    h(88) = K(22)*M(49);       % MP -- kDe --> M0
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% now turn the hazards into percentages with hw

    h_tot = sum(h);                 % Sum of hazard "bins"
    hw = h/h_tot;                   % Normalized "bins"
    
    % now pick a random number to step forward in time
    tt = -log(rand(1,1))/h_tot;
    t = t + tt;
    % pick a random number and use the weights to "select" an action
    r = rand(1,1);
    
    
    %% New
    if h_tot==0
        r = 10;
    end
    %% New
    
    % based on r, make a decision
    % Using M0
    if 0 <= r && r <= hw(1)
        M(1) = M(1) - 1;           
        X(1) = X(1) - 1;           
        M(2) = M(2) + 1;       
        reaction_count(1,runnum) = reaction_count(1,runnum) + 1;
        r_name = 1;
    elseif hw(1) < r && r <= sum(hw(1:2))
        M(2) = M(2) - 1;            
        M(1) = M(1) + 1;           
        X(1) = X(1) + 1;      
        reaction_count(2,runnum) = reaction_count(2,runnum) + 1;
        r_name = 2;
    elseif sum(hw(1:2)) < r && r <= sum(hw(1:3))
        M(2) = M(2) - 1;            
        M(3) = M(3) + 1; 
        reaction_count(3,runnum) = reaction_count(3,runnum) + 1;
        r_name = 3;
    elseif sum(hw(1:3)) < r && r <= sum(hw(1:4))
        M(3) = M(3) - 1;           
        M(4) = M(4) + 1; 
        reaction_count(4,runnum) = reaction_count(4,runnum) + 1;
        r_name = 4;
    elseif sum(hw(1:4)) < r && r <= sum(hw(1:5))
        M(4) = M(4) - 1;
        M(1) = M(1) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        reaction_count(5,runnum) = reaction_count(5,runnum) + 1;
        r_name = 5;
        
        % using M1
    elseif sum(hw(1:5)) < r && r <= sum(hw(1:6))
        M(6) = M(6) - 1;
        X(1) = X(1) - 1;
        M(7) = M(7) + 1;
        reaction_count(6,runnum) = reaction_count(6,runnum) + 1;
        r_name = 6;
    elseif sum(hw(1:6)) < r && r <= sum(hw(1:7))
        M(7) = M(7) - 1;
        M(6) = M(6) + 1;
        X(1) = X(1) + 1;
        reaction_count(7,runnum) = reaction_count(7,runnum) + 1;
        r_name = 7;
    elseif sum(hw(1:7)) < r && r <= sum(hw(1:8))
        M(7) = M(7) - 1;
        M(8) = M(8) + 1;
        reaction_count(8,runnum) = reaction_count(8,runnum) + 1;
        r_name = 8;
    elseif sum(hw(1:8)) < r && r <= sum(hw(1:9))
        M(8) = M(8) - 1;
        M(9) = M(9) + 1;
        reaction_count(9,runnum) = reaction_count(9,runnum) + 1;
        r_name = 9;
    elseif sum(hw(1:9)) < r && r <= sum(hw(1:10))
        M(9) = M(9) - 1;
        M(6) = M(6) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        reaction_count(10,runnum) = reaction_count(10,runnum) + 1;
        r_name = 10;
        
        % using M2
    elseif sum(hw(1:10)) < r && r <= sum(hw(1:11))
        M(13) = M(13) - 1;
        X(1)  = X(1) - 1;
        M(14) = M(14) + 1;
        reaction_count(11,runnum) = reaction_count(11,runnum) + 1;
        r_name = 11;
    elseif sum(hw(1:11)) < r && r <= sum(hw(1:12))
        M(14) = M(14) - 1;
        M(13) = M(13) + 1;
        X(1)  = X(1) + 1;
        reaction_count(12,runnum) = reaction_count(12,runnum) + 1;
        r_name = 12;
    elseif sum(hw(1:12)) < r && r <= sum(hw(1:13))
        M(14) = M(14) - 1;
        M(15) = M(15) + 1;
        reaction_count(13,runnum) = reaction_count(13,runnum) + 1;
        r_name = 13;
    elseif sum(hw(1:13)) < r && r <= sum(hw(1:14))
        M(15) = M(15) - 1;
        M(16) = M(16) + 1;
        reaction_count(14,runnum) = reaction_count(14,runnum) + 1;
        r_name = 14;
    elseif sum(hw(1:14)) < r && r <= sum(hw(1:15))
        M(16) = M(16) - 1;
        M(13) = M(13) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        reaction_count(15,runnum) = reaction_count(15,runnum) + 1;
        r_name = 15;
        
        % using M3
    elseif sum(hw(1:15)) < r && r <= sum(hw(1:16))
        M(20) = M(20) - 1;
        X(1) = X(1) - 1;
        M(21) = M(21) + 1;
        reaction_count(16,runnum) = reaction_count(16,runnum) + 1;
        r_name = 16;
    elseif sum(hw(1:16)) < r && r <= sum(hw(1:17))
        M(21) = M(21) - 1;
        M(20) = M(20) + 1;
        X(1) = X(1) + 1;
        reaction_count(17,runnum) = reaction_count(17,runnum) + 1;
        r_name = 17;
    elseif sum(hw(1:17)) < r && r <= sum(hw(1:18))
        M(21) = M(21) - 1;
        M(22) = M(22) + 1;
        reaction_count(18,runnum) = reaction_count(18,runnum) + 1;
        r_name = 18;
    elseif sum(hw(1:18)) < r && r <= sum(hw(1:19))
        M(22) = M(22) - 1;
        M(23) = M(23) + 1;
        reaction_count(19,runnum) = reaction_count(19,runnum) + 1;
        r_name = 19;
    elseif sum(hw(1:19)) < r && r <= sum(hw(1:20))
        M(23) = M(23) - 1;
        M(20) = M(20) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        reaction_count(20,runnum) = reaction_count(20,runnum) + 1;
        r_name = 20;
        
        % using M4
    elseif sum(hw(1:20)) < r && r <= sum(hw(1:21))
        M(27) = M(27) - 1;
        X(1) = X(1) - 1;
        M(28) = M(28) + 1;
        reaction_count(21,runnum) = reaction_count(21,runnum) + 1;
        r_name = 21;
    elseif sum(hw(1:21)) < r && r <= sum(hw(1:22))
        M(28) = M(28) - 1;
        M(27) = M(27) + 1;
        X(1) = X(1) + 1;
        reaction_count(22,runnum) = reaction_count(22,runnum) + 1;
        r_name = 22;
    elseif sum(hw(1:22)) < r && r <= sum(hw(1:23))
        M(28) = M(28) - 1;
        M(29) = M(29) + 1;
        reaction_count(23,runnum) = reaction_count(23,runnum) + 1;
        r_name = 23;
    elseif sum(hw(1:23)) < r && r <= sum(hw(1:24))
        M(29) = M(29) - 1;
        M(30) = M(30) + 1;
        reaction_count(24,runnum) = reaction_count(24,runnum) + 1;
        r_name = 24;
    elseif sum(hw(1:24)) < r && r <= sum(hw(1:25))
        M(30) = M(30) - 1;
        M(27) = M(27) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        reaction_count(25,runnum) = reaction_count(25,runnum) + 1;
        r_name = 25;
        
        % using M5
    elseif sum(hw(1:25)) < r && r <= sum(hw(1:26))
        M(34) = M(34) - 1;
        X(1) = X(1) - 1;
        M(35) = M(35) + 1;
        reaction_count(26,runnum) = reaction_count(26,runnum) + 1;
        r_name = 26;
    elseif sum(hw(1:26)) < r && r <= sum(hw(1:27))
        M(35) = M(35) - 1;
        M(34) = M(34) + 1;
        X(1) = X(1) + 1;
        reaction_count(27,runnum) = reaction_count(27,runnum) + 1;
        r_name = 27;
    elseif sum(hw(1:27)) < r && r <= sum(hw(1:28))
        M(35) = M(35) - 1;
        M(36) = M(36) + 1;
        reaction_count(28,runnum) = reaction_count(28,runnum) + 1;
        r_name = 28;
    elseif sum(hw(1:28)) < r && r <= sum(hw(1:29))
        M(36) = M(36) - 1;
        M(37) = M(37) + 1;
        reaction_count(29,runnum) = reaction_count(29,runnum) + 1;
        r_name = 29;
    elseif sum(hw(1:29)) < r && r <= sum(hw(1:30))
        M(37) = M(37) - 1;
        M(34) = M(34) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        reaction_count(30,runnum) = reaction_count(30,runnum) + 1;
        r_name = 30;
        
        % using M6
    elseif sum(hw(1:30)) < r && r <= sum(hw(1:31))
        M(41) = M(41) - 1;
        X(1) = X(1) - 1;
        M(42) = M(42) + 1;
        reaction_count(31,runnum) = reaction_count(31,runnum) + 1;
        r_name = 31;
    elseif sum(hw(1:31)) < r && r <= sum(hw(1:32))
        M(42) = M(42) - 1;
        M(41) = M(41) + 1;
        X(1) = X(1) + 1;
        reaction_count(32,runnum) = reaction_count(32,runnum) + 1;
        r_name = 32;
    elseif sum(hw(1:32)) < r && r <= sum(hw(1:33))
        M(42) = M(42) - 1;
        M(43) = M(43) + 1;
        reaction_count(33,runnum) = reaction_count(33,runnum) + 1;
        r_name = 33;
    elseif sum(hw(1:33)) < r && r <= sum(hw(1:34))
        M(43) = M(43) - 1;
        M(44) = M(44) + 1;
        reaction_count(34,runnum) = reaction_count(34,runnum) + 1;
        r_name = 34;
    elseif sum(hw(1:34)) < r && r <= sum(hw(1:35))
        M(44) = M(44) - 1;
        M(41) = M(41) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        reaction_count(35,runnum) = reaction_count(35,runnum) + 1;
        r_name = 35;
        
        % PLC and G-protein activation/deactivation
        
    elseif sum(hw(1:35)) < r && r <= sum(hw(1:36))
        X(2) = X(2) - 1;
        X(4) = X(4) - 1;
        X(5) = X(5) + 1;
        reaction_count(36,runnum) = reaction_count(36,runnum) + 1;
        r_name = 36;
    elseif sum(hw(1:36)) < r && r <= sum(hw(1:37))
        X(5) = X(5) - 1;
        X(9) = X(9) + 1;
        reaction_count(37,runnum) = reaction_count(37,runnum) + 1;
        r_name = 37;
    elseif sum(hw(1:73)) < r && r <= sum(hw(1:74))
        X(9) = X(9) - 1;
        X(4) = X(4) + 1;
        X(10) = X(10) + 1;
        reaction_count(74,runnum) = reaction_count(74,runnum) + 1;
        r_name = 74;
    elseif sum(hw(1:74)) < r && r <= sum(hw(1:75))
        X(10) = X(10) - 1;
        X(3) = X(3) - 1;
        X(1) = X(1) + 1;
        reaction_count(75,runnum) = reaction_count(75,runnum) + 1;
        r_name = 75;
        
        % Second Messenger Creation
        
    elseif sum(hw(1:37)) < r && r <= sum(hw(1:38))
        X(6) = X(6) + 1;
        reaction_count(38,runnum) = reaction_count(38,runnum) + 1;
        r_name = 38;
        
        % Channel Opening/Closing
        
    elseif sum(hw(1:38)) < r && r <= sum(hw(1:39))
        X(6) = X(6) - 1;
        X(7) = X(7) - 1;
        X(8) = X(8) + 1;
        reaction_count(39,runnum) = reaction_count(39,runnum) + 1;
        r_name = 39;
    elseif sum(hw(1:39)) < r && r <= sum(hw(1:40))
        X(8) = X(8) - 1;
        X(6) = X(6) + 1;
        X(7) = X(7) + 1;
        reaction_count(40,runnum) = reaction_count(40,runnum) + 1;
        r_name = 40;
        
        %Kinase Phosphorylation M0
        
    elseif sum(hw(1:40)) < r && r <= sum(hw(1:41))
        M(1) = M(1) - 1;
        M(5) = M(5) + 1;
        reaction_count(41,runnum) = reaction_count(41,runnum) + 1;
        r_name = 41;
    elseif sum(hw(1:41)) < r && r <= sum(hw(1:42))
        M(5) = M(5) - 1;
        M(1) = M(1) + 1;
        reaction_count(42,runnum) = reaction_count(42,runnum) + 1;
        r_name = 42;
    elseif sum(hw(1:42)) < r && r <= sum(hw(1:43))
        M(5) = M(5) - 1;
        M(10) = M(10) + 1;
        reaction_count(43,runnum) = reaction_count(43,runnum) + 1;
        r_name = 43;
        
        % M1
        
    elseif sum(hw(1:43)) < r && r <= sum(hw(1:44))
        M(6) = M(6) - 1;
        M(10) = M(10) + 1;
        reaction_count(44,runnum) = reaction_count(44,runnum) + 1;
        r_name = 44;
    elseif sum(hw(1:44)) < r && r <= sum(hw(1:45))
        M(10) = M(10) - 1;
        M(6) = M(6) + 1;
        reaction_count(45,runnum) = reaction_count(45,runnum) + 1;
        r_name = 45;
    elseif sum(hw(1:45)) < r && r <= sum(hw(1:46))
        M(10) = M(10) - 1;
        M(17) = M(17) + 1;
        reaction_count(46,runnum) = reaction_count(46,runnum) + 1;
        r_name = 46;
        
        %M2
        
    elseif sum(hw(1:46)) < r && r <= sum(hw(1:47))
        M(13) = M(13) - 1;
        M(17) = M(17) + 1;
        reaction_count(47,runnum) = reaction_count(47,runnum) + 1;
        r_name = 47;
    elseif sum(hw(1:47)) < r && r <= sum(hw(1:48))
        M(17) = M(17) - 1;
        M(13) = M(13) + 1;
        reaction_count(48,runnum) = reaction_count(48,runnum) + 1;
        r_name = 48;
    elseif sum(hw(1:48)) < r && r <= sum(hw(1:49))
        M(17) = M(17) - 1;
        M(24) = M(24) + 1;
        reaction_count(49,runnum) = reaction_count(49,runnum) + 1;
        r_name = 49;
        
        %M3
        
    elseif sum(hw(1:49)) < r && r <= sum(hw(1:50))
        M(20) = M(20) - 1;
        M(24) = M(24) + 1;
        reaction_count(50,runnum) = reaction_count(50,runnum) + 1;
        r_name = 50;
    elseif sum(hw(1:50)) < r && r <= sum(hw(1:51))
        M(24) = M(24) - 1;
        M(20) = M(20) + 1;
        reaction_count(51,runnum) = reaction_count(51,runnum) + 1;
        r_name = 51;
    elseif sum(hw(1:51)) < r && r <= sum(hw(1:52))
        M(24) = M(24) - 1;
        M(31) = M(31) + 1;
        reaction_count(52,runnum) = reaction_count(52,runnum) + 1;
        r_name = 52;
        
        %M4
        
    elseif sum(hw(1:52)) < r && r <= sum(hw(1:53))
        M(27) = M(27) - 1;
        M(31) = M(31) + 1;
        reaction_count(53,runnum) = reaction_count(53,runnum) + 1;
        r_name = 53;
    elseif sum(hw(1:53)) < r && r <= sum(hw(1:54))
        M(31) = M(31) - 1;
        M(27) = M(27) + 1;
        reaction_count(54,runnum) = reaction_count(54,runnum) + 1;
        r_name = 54;
    elseif sum(hw(1:54)) < r && r <= sum(hw(1:55))
        M(31) = M(31) - 1;
        M(38) = M(38) + 1;
        reaction_count(55,runnum) = reaction_count(55,runnum) + 1;
        r_name = 55;
        
        %M5
        
    elseif sum(hw(1:55)) < r && r <= sum(hw(1:56))
        M(34) = M(34) - 1;
        M(38) = M(38) + 1;
        reaction_count(56,runnum) = reaction_count(56,runnum) + 1;
        r_name = 56;
    elseif sum(hw(1:56)) < r && r <= sum(hw(1:57))
        M(38) = M(38) - 1;
        M(34) = M(34) + 1;
        reaction_count(57,runnum) = reaction_count(57,runnum) + 1;
        r_name = 57;
    elseif sum(hw(1:57)) < r && r <= sum(hw(1:58))
        M(38) = M(38) - 1;
        M(45) = M(45) + 1;
        reaction_count(58,runnum) = reaction_count(58,runnum) + 1;
        r_name = 58;
        
        %M6
        
    elseif sum(hw(1:58)) < r && r <= sum(hw(1:59))
        M(41) = M(41) - 1;
        M(45) = M(45) + 1;
        reaction_count(59,runnum) = reaction_count(59,runnum) + 1;
        r_name = 59;
    elseif sum(hw(1:59)) < r && r <= sum(hw(1:60))
        M(45) = M(45) - 1;
        M(41) = M(41) + 1;
        reaction_count(60,runnum) = reaction_count(60,runnum) + 1;
        r_name = 60;
        
        % Arrestin Binding
        % M1      
    elseif sum(hw(1:60)) < r && r <= sum(hw(1:61))
        M(6) = M(6) - 1;
        X(11) = X(11) - 1;
        M(11) = M(11) + 1;
        reaction_count(61,runnum) = reaction_count(61,runnum) + 1;
        r_name = 61;
    elseif sum(hw(1:61)) < r && r <= sum(hw(1:62))
        M(6) = M(6) - 1;
        X(12) = X(12) - 1;
        M(12) = M(12) + 1;
        reaction_count(62,runnum) = reaction_count(62,runnum) + 1;
        r_name = 62;
        
        %M2
    elseif sum(hw(1:62)) < r && r <= sum(hw(1:63))
        M(13) = M(13) - 1;
        X(11) = X(11) - 1;
        M(18) = M(18) + 1;
        reaction_count(63,runnum) = reaction_count(63,runnum) + 1;
        r_name = 63;
    elseif sum(hw(1:63)) < r && r <= sum(hw(1:64))
        M(13) = M(13) - 1;
        X(12) = X(12) - 1;
        M(19) = M(19) + 1;
        reaction_count(64,runnum) = reaction_count(64,runnum) + 1;
        r_name = 64;
        
        %M3       
    elseif sum(hw(1:64)) < r && r <= sum(hw(1:65))
        M(20) = M(20) - 1;
        X(11) = X(11) - 1;
        M(25) = M(25) + 1;
        reaction_count(65,runnum) = reaction_count(65,runnum) + 1;
        r_name = 65;
    elseif sum(hw(1:65)) < r && r <= sum(hw(1:66))
        M(20) = M(20) - 1;
        X(12) = X(12) - 1;
        M(26) = M(26) + 1;
        reaction_count(66,runnum) = reaction_count(66,runnum) + 1;
        r_name = 66;
        
        %M4      
    elseif sum(hw(1:66)) < r && r <= sum(hw(1:67))
        M(27) = M(27) - 1;
        X(11) = X(11) - 1;
        M(32) = M(32) + 1;
        reaction_count(67,runnum) = reaction_count(67,runnum) + 1;
        r_name = 67;
    elseif sum(hw(1:67)) < r && r <= sum(hw(1:68))
        M(27) = M(27) - 1;
        X(12) = X(12) - 1;
        M(33) = M(33) + 1;
        reaction_count(68,runnum) = reaction_count(68,runnum) + 1;
        r_name = 68;
        
        %M5
    elseif sum(hw(1:68)) < r && r <= sum(hw(1:69))
        M(34) = M(34) - 1;
        X(11) = X(11) - 1;
        M(39) = M(39) + 1;
        reaction_count(69,runnum) = reaction_count(69,runnum) + 1;
        r_name = 69;
    elseif sum(hw(1:69)) < r && r <= sum(hw(1:70))
        M(34) = M(34) - 1;
        X(12) = X(12) - 1;
        M(40) = M(40) + 1;
        reaction_count(70,runnum) = reaction_count(70,runnum) + 1;
        r_name = 70;
        
        %M6
    elseif sum(hw(1:70)) < r && r <= sum(hw(1:71))
        M(41) = M(41) - 1;
        X(11) = X(11) - 1;
        M(46) = M(46) + 1;
        reaction_count(71,runnum) = reaction_count(71,runnum) + 1;
        r_name = 71;
    elseif sum(hw(1:71)) < r && r <= sum(hw(1:72))
        M(41) = M(41) - 1;
        X(12) = X(12) - 1;
        M(47) = M(47) + 1;
        reaction_count(72,runnum) = reaction_count(72,runnum) + 1;
        r_name = 72;
        
        % SecM degradation
    elseif sum(hw(1:72)) < r && r <= sum(hw(1:73))
        X(6) = X(6) - 1;
        reaction_count(73,runnum) = reaction_count(73,runnum) + 1;
        r_name = 73;
    
        % Arresin unbinding 
        
       
    elseif sum(hw(1:75)) < r && r <= sum(hw(1:76))
         M(11) = M(11) - 1;
         X(11) = X(11) + 1;
         M(49) = M(49) + 1;
         reaction_count(76,runnum) = reaction_count(76,runnum) + 1;
         r_name = 76;
       
    elseif sum(hw(1:76)) < r && r <= sum(hw(1:77))
         M(18) = M(18) - 1;
         X(11) = X(11) + 1;
         M(49) = M(49) + 1;
         reaction_count(77,runnum) = reaction_count(77,runnum) + 1;
         r_name = 77;
       
    elseif sum(hw(1:77)) < r && r <= sum(hw(1:78))     
         M(25) = M(25) - 1;
         X(11) = X(11) + 1;
         M(49) = M(49) + 1;
         reaction_count(78,runnum) = reaction_count(78,runnum) + 1;
         r_name = 78;
        
    elseif sum(hw(1:78)) < r && r <= sum(hw(1:79))
         M(32) = M(32) - 1;
         X(11) = X(11) + 1;
         M(49) = M(49) + 1;
         reaction_count(79,runnum) = reaction_count(79,runnum) + 1;
         r_name = 79;
        
    elseif sum(hw(1:79)) < r && r <= sum(hw(1:80))     
         M(39) = M(39) - 1;
         X(11) = X(11) + 1;
         M(49) = M(49) + 1;
         reaction_count(80,runnum) = reaction_count(80,runnum) + 1;
         r_name = 80;
        
    elseif sum(hw(1:80)) < r && r <= sum(hw(1:81))      
         M(46) = M(46) - 1;
         X(11) = X(11) + 1;
         M(49) = M(49) + 1;
         reaction_count(81,runnum) = reaction_count(81,runnum) + 1;
         r_name = 81;
        
    elseif sum(hw(1:81)) < r && r <= sum(hw(1:82))
         M(12) = M(12) - 1;
         X(12) = X(12) + 1;
         M(49) = M(49) + 1;
         reaction_count(82,runnum) = reaction_count(82,runnum) + 1;
         r_name = 82;
         
    elseif sum(hw(1:82)) < r && r <= sum(hw(1:83))
         M(19) = M(19) - 1;
         X(12) = X(12) + 1;
         M(49) = M(49) + 1;
         reaction_count(83,runnum) = reaction_count(83,runnum) + 1;
         r_name = 83;
  
    elseif sum(hw(1:83)) < r && r <= sum(hw(1:84))      
         M(26) = M(26) - 1;
         X(12) = X(12) + 1;
         M(49) = M(49) + 1;
         reaction_count(84,runnum) = reaction_count(84,runnum) + 1;
         r_name = 84;

    elseif sum(hw(1:84)) < r && r <= sum(hw(1:85))
         M(33) = M(33) - 1;
         X(12) = X(12) + 1;
         M(49) = M(49) + 1;
         reaction_count(85,runnum) = reaction_count(85,runnum) + 1;
         r_name = 85;
 
    elseif sum(hw(1:85)) < r && r <= sum(hw(1:86))    
         M(40) = M(40) - 1;
         X(12) = X(12) + 1;
         M(49) = M(49) + 1;
         reaction_count(86,runnum) = reaction_count(86,runnum) + 1;
         r_name = 86;
 
    elseif sum(hw(1:86)) < r && r <= sum(hw(1:87))
         M(47) = M(47) - 1;
         X(12) = X(12) + 1;
         M(49) = M(49) + 1;
         reaction_count(87,runnum) = reaction_count(87,runnum) + 1;
         r_name = 87;
 
    elseif sum(hw(1:87)) < r && r <= sum(hw(1:88))
         M(48) = M(48) + 1;
         M(49) = M(49) - 1;
         reaction_count(88,runnum) = reaction_count(88,runnum) + 1;
         r_name = 88;
         
  
         else
        M(1) = M(1);   % If the # of melanopsin cells is zero
   
    end
    
    %% beginning, graph with Time
    Time = [0:time_step:tmax]; 
    %% beginning, graph with Time
    
    %% store time, molecule numbers in every 'time_step' sec
    %% if h_tot==0, move time t to tflash
    if h_tot==0
        t = min(tflash,tmax);
    end
    %% if h_tot==0, move time t to tflash
    
     time_index = floor(t/time_step) + 1; % an index for the current time
    
    %% if time_index becomes larger than size of tstore 
    %% (final time can be larger than tmax), 
    %% make time_index equal to the the final index 
    if time_index > floor(tmax/time_step)+1;
        time_index = floor(tmax/time_step)+1;
    end
    %% 
    reaction_sequence(counter,runnum) = r_name;
    
    %% store time, molecule numbers in every 'time_step' sec
    
    if time_index > prev_t_index
        for j = (prev_t_index+1):time_index
            tstore(j,1,runnum) = t;
            Xstore(j,:,runnum) = X;
            Mstore(j,:,runnum) = M;
            ttstore(j,1,runnum) = tt;
        end
    end
    prev_t_index = time_index; 
    %% store time, molecule numbers in every 'time_step' sec

%    whatha=whatha+1 
% %     if t>19.89
% %         keyboard;
% %     end
% % %    pause(0.2)
% %     plot_matrx(tic,:) = [t X M];
% % figure(1)    
% %     % if the new t < tmax, repeat.  Else plot the results.
% %     subplot(2,1,1)
% %     plot(M)
% %     pause(0.1)
% %     subplot(2,1,2)
% %     plot(X)
% %    axis([0 14 0 100])
% %    pause(0.1)
% %    display(X(1))
   %%display(X(5))
   %%pause(1)

end % end of one realization
% simstuff = {tstore,Mstore,Xstore};
% save(sprintf('run%d',runnum),'simstuff')
save('results.mat','tstore','Time','Mstore','Xstore','reaction_count','seed_name','reaction_sequence');

end % end of all realizations

%% compute mean and standard deviation
% Mx = mean(Xstore,3);
% Sx = std(Xstore,0,3);
% Mm = mean(Mstore,3);
% Sm = std(Mstore,0,3);

%% opchan = Xstore(:,8,:)./(Xstore(:,7,:)+Xstore(:,8,:));
%% Instead of using mean and std of opchan, we can use those of 
%% Xstore(:,8,:) since Xstore(:,7,:)+Xstore(:,8,:)is constant
Mopchan = mean(Xstore(:,8,:),3);
% '3' means averaging with respect to the 3rd dimension.
Sopchan = std(Xstore(:,8,:),0,3);
%% Scaling everything by max(Mopchan) (appropriate?)
Mratio = Mopchan./max(Mopchan);
Sratio = Sopchan./max(Mopchan);
%% compute mean and standard deviation

save('results.mat','Time','tstore','Mstore','Xstore','Mopchan','Sopchan','Mratio','Sratio','reaction_count','seed_name','reaction_sequence');


%% plotting average
% % figure(1)
% % plot(Time,Mm(:,1)+Mm(:,6)+Mm(:,13)+Mm(:,20)+Mm(:,27)+Mm(:,34)+Mm(:,41));
% % %axis([0 tmax -1 100]);
% % xlabel('time (/sec)'); ylabel('# of cells');
LowPassFilter
Exp_data = [newdata(1100:end,1)-1.1 newdata(1100:end,2)]

% Exp_data =dlmread('ModifiedData.csv'); % calcium imaging
% Column1: Time; 
% Column2: Mel mean; Column3: P-null mean; Column4: P2St/Mel mean; 
% Column5: P2St/P-null mean; Column6: P1St/Mel mean;
% Column7: Mel std; Column8: P-null std; Column9: P2St/Mel std; 
% Column10: P2St/P-null std; Column11: P1St/Mel std;
Data_Time = Exp_data(:,1);
Data_Mean = Exp_data(:,2);
% Data_Std = Exp_data(:,7:11);

figure(2)
hold on
% errorbar(Data_Time,Data_Mean(:,1),Data_Std(:,1),'kx') 
plot(Data_Time,Data_Mean(:,1),'r-','LineWidth',4)
set(gcf,'color','w');
grid on
axis([0 5 -0.05 1.4])
box on

% errorbar(Time,Mratio,Sratio,'bx') 
%plot(Time,Mratio+Sratio,'b:','LineWidth',1)
%plot(Time,Mratio-Sratio,'b:','LineWidth',1)
%% Comment the above line if we hide standard deviations
plot(Time,Mratio,'g-','LineWidth',4)
% a ratio of the number of open channels out of the total number of
% channels
hold off
%% plotting average


%% plot mean of reaction_count
clear hh
hh = mean(reaction_count,2);
GproteinActivation = [hh(1) hh(6) hh(11) hh(16) hh(2) hh(7) hh(12) hh(17) hh(3) hh(8) ...
    hh(13) hh(18) hh(4) hh(9) hh(14) hh(19) hh(5) hh(10) hh(15) hh(20)];
KinasePhosphorylation = [hh(41) hh(44) hh(47) hh(50) hh(42) hh(45) hh(48) hh(51) hh(43) hh(46) hh(49)];
%% check KinasePhosphorylation
ArrestinBindingUnbinding = [hh(61) hh(63) hh(65) hh(62) hh(64) hh(66) hh(76) hh(77) hh(78) hh(82) hh(83) hh(84) h(88)];
PLCSecMChannel = [hh(36) hh(37) hh(74) hh(75) hh(38) hh(73) hh(39) hh(40)];   

figure(3)
subplot(4,1,1)
hold on
bar([1:1:size(GproteinActivation,2)],GproteinActivation)
set(gca,'XTick',[1:1:size(GproteinActivation,2)])
set(gca,'XTickLabel',{'kG1y(0)';'kG1y(1)';'kG1y(2)';'kG1y(3)';'kG2(0)';'kG2(1)';'kG2(2)';'kG2(3)';...
    'kG3(0)';'kG3(1)';'kG3(2)';'kG3(3)';'kG4(0)';'kG4(1)';'kG4(2)';'kG4(3)';'kG5(0)';'kG5(1)';'kG5(2)';'kG5(3)'})
set(gcf,'color','w');
grid on
box on
hold off

subplot(4,1,2)
hold on
bar([1:1:size(KinasePhosphorylation,2)],KinasePhosphorylation')
set(gca,'XTick',[1:1:size(KinasePhosphorylation,2)])
set(gca,'XTickLabel',{'kk1(0)';'kk1(1)';'kk1(2)';'kk1(3)';'kk2(0)';'kk2(1)';'kk2(2)';'kk2(3)';'kk3(0)';'kk3(1)';'kk3(2)'})
set(gcf,'color','w');
grid on
box on
hold off

subplot(4,1,3)
hold on
bar([1:1:size(ArrestinBindingUnbinding,2)],ArrestinBindingUnbinding')
set(gca,'XTick',[1:1:size(ArrestinBindingUnbinding,2)])
set(gca,'XTickLabel',{'kB1w(1)';'kB1w(2)';'kB1w(3)';'kB2w(1)';'kB2w(2)';'kB2w(3)';'kUB1w(1)';'kUB1w(2)';'kUB1w(3)';'kUB2w(1)';'kUB2w(2)';'kUB2w(3)';'kDe'})
set(gcf,'color','w');
grid on
box on
hold off  

subplot(4,1,4)
hold on
bar([1:1:size(PLCSecMChannel,2)],PLCSecMChannel')
set(gca,'XTick',[1:1:size(PLCSecMChannel,2)])
set(gca,'XTickLabel',{'kP';'kI1';'kI2';'kI3';'kS';'delta';'kO';'kC'})
set(gcf,'color','w');
grid on
box on
hold off
%% plot mean of reaction_count
%% check the label one more time and send to Christy

toc
clear variables