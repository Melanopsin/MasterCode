% Melanopsin_Sensitivity(totrun,dataset) will run the Gillespie algorithm code
% for one run and using the initial values 
% and stochastic rate constants for the specified dataset (dataset). Later once
% the optimization code is finished the rate constants used will be from
% the rates written to file by that routine. 

function [Time,open_channels] = Melanopsin_Sensitivity(tmax,flashint,LHSmatrix,x)
%% store time, molecule numbers in every 'time_step' sec for all runs
time_step=0.25; % maximum time step, before we set as 0.25
maxcounter=10000000;

%% You record t, M, X in every time interval equal to time_step
Xstore = zeros(floor(tmax/time_step)+1,12); 
Mstore = zeros(floor(tmax/time_step)+1,49); 
ttstore = zeros(floor(tmax/time_step)+1,1);

no_rxns = 88;                   % (<--- this is hardcode) number of reactions (total)
%% store time, molecule numbers in every 'time_step' sec for all runs
slowmedown=1.00;
tic;

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


%% set parameters and initial conditions in each run
Melanopsin_Parameter_settings_LHS;

GTP = LHSmatrix(x,1);
PIP2 = LHSmatrix(x,2);
Ki = LHSmatrix(x,3);
ATP = LHSmatrix(x,4);
kmax = LHSmatrix(x,5);
KM = LHSmatrix(x,6);

kG1 = LHSmatrix(x,7); 
kG2 = LHSmatrix(x,8); 
kG3 = LHSmatrix(x,9); 
kG4 = LHSmatrix(x,10);
kG5 = LHSmatrix(x,11);

kP = LHSmatrix(x,12);
kI1 = LHSmatrix(x,13);
kI2 = LHSmatrix(x,14);
kI3 = LHSmatrix(x,15);
kS = LHSmatrix(x,16);
kO = LHSmatrix(x,17);
kC = LHSmatrix(x,18); 

kk1 = LHSmatrix(x,19); 
kk2 = LHSmatrix(x,20); 
kk3 = LHSmatrix(x,21);
kB1 = LHSmatrix(x,22);
kB2 = LHSmatrix(x,23);
kUB1 = LHSmatrix(x,24);
kUB2 = LHSmatrix(x,25);
kDe = LHSmatrix(x,26);


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
h = zeros(no_rxns,1);  % initialize the hazard vector
tstore(1,1) = t;
Xstore(1,:) = X;
Mstore(1,:) = M;
ttstore(1,1) = 0;
prev_t_index = 1; % it stores the previous time intex
%% set initial values


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
    h(21) = 0*Y(4)*K(1)*M(27)*X(1);
    h(26) = 0*Y(5)*K(1)*M(34)*X(1);
    h(31) = 0*Y(6)*K(1)*M(41)*X(1);

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
    elseif hw(1) < r && r <= sum(hw(1:2))
        M(2) = M(2) - 1;            
        M(1) = M(1) + 1;           
        X(1) = X(1) + 1;      
    elseif sum(hw(1:2)) < r && r <= sum(hw(1:3))
        M(2) = M(2) - 1;            
        M(3) = M(3) + 1; 
    elseif sum(hw(1:3)) < r && r <= sum(hw(1:4))
        M(3) = M(3) - 1;           
        M(4) = M(4) + 1; 
    elseif sum(hw(1:4)) < r && r <= sum(hw(1:5))
        M(4) = M(4) - 1;
        M(1) = M(1) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        
        % using M1
    elseif sum(hw(1:5)) < r && r <= sum(hw(1:6))
        M(6) = M(6) - 1;
        X(1) = X(1) - 1;
        M(7) = M(7) + 1;
    elseif sum(hw(1:6)) < r && r <= sum(hw(1:7))
        M(7) = M(7) - 1;
        M(6) = M(6) + 1;
        X(1) = X(1) + 1;
    elseif sum(hw(1:7)) < r && r <= sum(hw(1:8))
        M(7) = M(7) - 1;
        M(8) = M(8) + 1;
    elseif sum(hw(1:8)) < r && r <= sum(hw(1:9))
        M(8) = M(8) - 1;
        M(9) = M(9) + 1;
    elseif sum(hw(1:9)) < r && r <= sum(hw(1:10))
        M(9) = M(9) - 1;
        M(6) = M(6) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        
        % using M2
    elseif sum(hw(1:10)) < r && r <= sum(hw(1:11))
        M(13) = M(13) - 1;
        X(1)  = X(1) - 1;
        M(14) = M(14) + 1;
    elseif sum(hw(1:11)) < r && r <= sum(hw(1:12))
        M(14) = M(14) - 1;
        M(13) = M(13) + 1;
        X(1)  = X(1) + 1;
    elseif sum(hw(1:12)) < r && r <= sum(hw(1:13))
        M(14) = M(14) - 1;
        M(15) = M(15) + 1;
    elseif sum(hw(1:13)) < r && r <= sum(hw(1:14))
        M(15) = M(15) - 1;
        M(16) = M(16) + 1;
    elseif sum(hw(1:14)) < r && r <= sum(hw(1:15))
        M(16) = M(16) - 1;
        M(13) = M(13) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        
        % using M3
    elseif sum(hw(1:15)) < r && r <= sum(hw(1:16))
        M(20) = M(20) - 1;
        X(1) = X(1) - 1;
        M(21) = M(21) + 1;
    elseif sum(hw(1:16)) < r && r <= sum(hw(1:17))
        M(21) = M(21) - 1;
        M(20) = M(20) + 1;
        X(1) = X(1) + 1;
    elseif sum(hw(1:17)) < r && r <= sum(hw(1:18))
        M(21) = M(21) - 1;
        M(22) = M(22) + 1;
    elseif sum(hw(1:18)) < r && r <= sum(hw(1:19))
        M(22) = M(22) - 1;
        M(23) = M(23) + 1;
    elseif sum(hw(1:19)) < r && r <= sum(hw(1:20))
        M(23) = M(23) - 1;
        M(20) = M(20) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        
        % using M4
    elseif sum(hw(1:20)) < r && r <= sum(hw(1:21))
        M(27) = M(27) - 1;
        X(1) = X(1) - 1;
        M(28) = M(28) + 1;
    elseif sum(hw(1:21)) < r && r <= sum(hw(1:22))
        M(28) = M(28) - 1;
        M(27) = M(27) + 1;
        X(1) = X(1) + 1;
    elseif sum(hw(1:22)) < r && r <= sum(hw(1:23))
        M(28) = M(28) - 1;
        M(29) = M(29) + 1;
    elseif sum(hw(1:23)) < r && r <= sum(hw(1:24))
        M(29) = M(29) - 1;
        M(30) = M(30) + 1;
    elseif sum(hw(1:24)) < r && r <= sum(hw(1:25))
        M(30) = M(30) - 1;
        M(27) = M(27) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        
        % using M5
    elseif sum(hw(1:25)) < r && r <= sum(hw(1:26))
        M(34) = M(34) - 1;
        X(1) = X(1) - 1;
        M(35) = M(35) + 1;
    elseif sum(hw(1:26)) < r && r <= sum(hw(1:27))
        M(35) = M(35) - 1;
        M(34) = M(34) + 1;
        X(1) = X(1) + 1;
    elseif sum(hw(1:27)) < r && r <= sum(hw(1:28))
        M(35) = M(35) - 1;
        M(36) = M(36) + 1;
    elseif sum(hw(1:28)) < r && r <= sum(hw(1:29))
        M(36) = M(36) - 1;
        M(37) = M(37) + 1;
    elseif sum(hw(1:29)) < r && r <= sum(hw(1:30))
        M(37) = M(37) - 1;
        M(34) = M(34) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        
        % using M6
    elseif sum(hw(1:30)) < r && r <= sum(hw(1:31))
        M(41) = M(41) - 1;
        X(1) = X(1) - 1;
        M(42) = M(42) + 1;
    elseif sum(hw(1:31)) < r && r <= sum(hw(1:32))
        M(42) = M(42) - 1;
        M(41) = M(41) + 1;
        X(1) = X(1) + 1;
    elseif sum(hw(1:32)) < r && r <= sum(hw(1:33))
        M(42) = M(42) - 1;
        M(43) = M(43) + 1;
    elseif sum(hw(1:33)) < r && r <= sum(hw(1:34))
        M(43) = M(43) - 1;
        M(44) = M(44) + 1;
    elseif sum(hw(1:34)) < r && r <= sum(hw(1:35))
        M(44) = M(44) - 1;
        M(41) = M(41) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        
        % PLC and G-protein activation/deactivation
        
    elseif sum(hw(1:35)) < r && r <= sum(hw(1:36))
        X(2) = X(2) - 1;
        X(4) = X(4) - 1;
        X(5) = X(5) + 1;
    elseif sum(hw(1:36)) < r && r <= sum(hw(1:37))
        X(5) = X(5) - 1;
        X(9) = X(9) + 1;
    elseif sum(hw(1:73)) < r && r <= sum(hw(1:74))
        X(9) = X(9) - 1;
        X(4) = X(4) + 1;
        X(10) = X(10) + 1;
    elseif sum(hw(1:74)) < r && r <= sum(hw(1:75))
        X(10) = X(10) - 1;
        X(3) = X(3) - 1;
        X(1) = X(1) + 1;
        
        % Second Messenger Creation
        
    elseif sum(hw(1:37)) < r && r <= sum(hw(1:38))
        X(6) = X(6) + 1;
        
        % Channel Opening/Closing
        
    elseif sum(hw(1:38)) < r && r <= sum(hw(1:39))
        X(6) = X(6) - 1;
        X(7) = X(7) - 1;
        X(8) = X(8) + 1;
    elseif sum(hw(1:39)) < r && r <= sum(hw(1:40))
        X(8) = X(8) - 1;
        X(6) = X(6) + 1;
        X(7) = X(7) + 1;
        
        %Kinase Phosphorylation M0
        
    elseif sum(hw(1:40)) < r && r <= sum(hw(1:41))
        M(1) = M(1) - 1;
        M(5) = M(5) + 1;
    elseif sum(hw(1:41)) < r && r <= sum(hw(1:42))
        M(5) = M(5) - 1;
        M(1) = M(1) + 1;
    elseif sum(hw(1:42)) < r && r <= sum(hw(1:43))
        M(5) = M(5) - 1;
        M(10) = M(10) + 1;
        
        % M1
        
    elseif sum(hw(1:43)) < r && r <= sum(hw(1:44))
        M(6) = M(6) - 1;
        M(10) = M(10) + 1;
    elseif sum(hw(1:44)) < r && r <= sum(hw(1:45))
        M(10) = M(10) - 1;
        M(6) = M(6) + 1;
    elseif sum(hw(1:45)) < r && r <= sum(hw(1:46))
        M(10) = M(10) - 1;
        M(17) = M(17) + 1;
        
        %M2
        
    elseif sum(hw(1:46)) < r && r <= sum(hw(1:47))
        M(13) = M(13) - 1;
        M(17) = M(17) + 1;
    elseif sum(hw(1:47)) < r && r <= sum(hw(1:48))
        M(17) = M(17) - 1;
        M(13) = M(13) + 1;
    elseif sum(hw(1:48)) < r && r <= sum(hw(1:49))
        M(17) = M(17) - 1;
        M(24) = M(24) + 1;
        
        %M3
        
    elseif sum(hw(1:49)) < r && r <= sum(hw(1:50))
        M(20) = M(20) - 1;
        M(24) = M(24) + 1;
    elseif sum(hw(1:50)) < r && r <= sum(hw(1:51))
        M(24) = M(24) - 1;
        M(20) = M(20) + 1;
    elseif sum(hw(1:51)) < r && r <= sum(hw(1:52))
        M(24) = M(24) - 1;
        M(31) = M(31) + 1;
        
        %M4
        
    elseif sum(hw(1:52)) < r && r <= sum(hw(1:53))
        M(27) = M(27) - 1;
        M(31) = M(31) + 1;
    elseif sum(hw(1:53)) < r && r <= sum(hw(1:54))
        M(31) = M(31) - 1;
        M(27) = M(27) + 1;
    elseif sum(hw(1:54)) < r && r <= sum(hw(1:55))
        M(31) = M(31) - 1;
        M(38) = M(38) + 1;
        
        %M5
        
    elseif sum(hw(1:55)) < r && r <= sum(hw(1:56))
        M(34) = M(34) - 1;
        M(38) = M(38) + 1;
    elseif sum(hw(1:56)) < r && r <= sum(hw(1:57))
        M(38) = M(38) - 1;
        M(34) = M(34) + 1;
    elseif sum(hw(1:57)) < r && r <= sum(hw(1:58))
        M(38) = M(38) - 1;
        M(45) = M(45) + 1;
        
        %M6
        
    elseif sum(hw(1:58)) < r && r <= sum(hw(1:59))
        M(41) = M(41) - 1;
        M(45) = M(45) + 1;
    elseif sum(hw(1:59)) < r && r <= sum(hw(1:60))
        M(45) = M(45) - 1;
        M(41) = M(41) + 1;
        
        % Arrestin Binding
        % M1      
    elseif sum(hw(1:60)) < r && r <= sum(hw(1:61))
        M(6) = M(6) - 1;
        X(11) = X(11) - 1;
        M(11) = M(11) + 1;
    elseif sum(hw(1:61)) < r && r <= sum(hw(1:62))
        M(6) = M(6) - 1;
        X(12) = X(12) - 1;
        M(12) = M(12) + 1;
        
        %M2
    elseif sum(hw(1:62)) < r && r <= sum(hw(1:63))
        M(13) = M(13) - 1;
        X(11) = X(11) - 1;
        M(18) = M(18) + 1;
    elseif sum(hw(1:63)) < r && r <= sum(hw(1:64))
        M(13) = M(13) - 1;
        X(12) = X(12) - 1;
        M(19) = M(19) + 1;
        
        %M3       
    elseif sum(hw(1:64)) < r && r <= sum(hw(1:65))
        M(20) = M(20) - 1;
        X(11) = X(11) - 1;
        M(25) = M(25) + 1;
    elseif sum(hw(1:65)) < r && r <= sum(hw(1:66))
        M(20) = M(20) - 1;
        X(12) = X(12) - 1;
        M(26) = M(26) + 1;
        
        %M4      
    elseif sum(hw(1:66)) < r && r <= sum(hw(1:67))
        M(27) = M(27) - 1;
        X(11) = X(11) - 1;
        M(32) = M(32) + 1;
    elseif sum(hw(1:67)) < r && r <= sum(hw(1:68))
        M(27) = M(27) - 1;
        X(12) = X(12) - 1;
        M(33) = M(33) + 1;
        
        %M5
    elseif sum(hw(1:68)) < r && r <= sum(hw(1:69))
        M(34) = M(34) - 1;
        X(11) = X(11) - 1;
        M(39) = M(39) + 1;
    elseif sum(hw(1:69)) < r && r <= sum(hw(1:70))
        M(34) = M(34) - 1;
        X(12) = X(12) - 1;
        M(40) = M(40) + 1;
        
        %M6
    elseif sum(hw(1:70)) < r && r <= sum(hw(1:71))
        M(41) = M(41) - 1;
        X(11) = X(11) - 1;
        M(46) = M(46) + 1;
    elseif sum(hw(1:71)) < r && r <= sum(hw(1:72))
        M(41) = M(41) - 1;
        X(12) = X(12) - 1;
        M(47) = M(47) + 1;
        
        % SecM degradation
    elseif sum(hw(1:72)) < r && r <= sum(hw(1:73))
        X(6) = X(6) - 1;
    
        % Arresin unbinding 
        
       
    elseif sum(hw(1:75)) < r && r <= sum(hw(1:76))
         M(11) = M(11) - 1;
         X(11) = X(11) + 1;
         M(49) = M(49) + 1;
       
    elseif sum(hw(1:76)) < r && r <= sum(hw(1:77))
         M(18) = M(18) - 1;
         X(11) = X(11) + 1;
         M(49) = M(49) + 1;
       
    elseif sum(hw(1:77)) < r && r <= sum(hw(1:78))     
         M(25) = M(25) - 1;
         X(11) = X(11) + 1;
         M(49) = M(49) + 1;
        
    elseif sum(hw(1:78)) < r && r <= sum(hw(1:79))
         M(32) = M(32) - 1;
         X(11) = X(11) + 1;
         M(49) = M(49) + 1;
        
    elseif sum(hw(1:79)) < r && r <= sum(hw(1:80))     
         M(39) = M(39) - 1;
         X(11) = X(11) + 1;
         M(49) = M(49) + 1;
        
    elseif sum(hw(1:80)) < r && r <= sum(hw(1:81))      
         M(46) = M(46) - 1;
         X(11) = X(11) + 1;
         M(49) = M(49) + 1;
        
    elseif sum(hw(1:81)) < r && r <= sum(hw(1:82))
         M(12) = M(12) - 1;
         X(12) = X(12) + 1;
         M(49) = M(49) + 1;
         
    elseif sum(hw(1:82)) < r && r <= sum(hw(1:83))
         M(19) = M(19) - 1;
         X(12) = X(12) + 1;
         M(49) = M(49) + 1;
  
    elseif sum(hw(1:83)) < r && r <= sum(hw(1:84))      
         M(26) = M(26) - 1;
         X(12) = X(12) + 1;
         M(49) = M(49) + 1;

    elseif sum(hw(1:84)) < r && r <= sum(hw(1:85))
         M(33) = M(33) - 1;
         X(12) = X(12) + 1;
         M(49) = M(49) + 1;
 
    elseif sum(hw(1:85)) < r && r <= sum(hw(1:86))    
         M(40) = M(40) - 1;
         X(12) = X(12) + 1;
         M(49) = M(49) + 1;
 
    elseif sum(hw(1:86)) < r && r <= sum(hw(1:87))
         M(47) = M(47) - 1;
         X(12) = X(12) + 1;
         M(49) = M(49) + 1;
 
    elseif sum(hw(1:87)) < r && r <= sum(hw(1:88))
         M(48) = M(48) + 1;
         M(49) = M(49) - 1;
         
  
         else
        M(1) = M(1);   % If the # of melanopsin cells is zero
   
    end
    
    %% beginning, graph with Time
    Time = (0:time_step:tmax)'; 
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
    
    %% store time, molecule numbers in every 'time_step' sec
    
    if time_index > prev_t_index
        for j = (prev_t_index+1):time_index
            tstore(j,1) = t;
            Xstore(j,:) = X;
            Mstore(j,:) = M;
            ttstore(j,1) = tt;
        end
    end
    prev_t_index = time_index; 
    %% store time, molecule numbers in every 'time_step' sec

end % end of one realization

open_channels = Xstore(:,8);
