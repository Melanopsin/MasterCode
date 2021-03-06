% Melanopsin2_count will run a single run of the same code found in melanopsin.m
% as Melanopsin2 does. 'reaction_count' records each reaction's occurence so that 
% we can check whether every reaction happens.

% function [tout Mout Xout max_chan maxchan_time] = Melanopsin2_count(dataset)
function Melanopsin2_count(dataset)

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
% Mn* + ArrB1 -- kK4 w(n) --> Mn*.ArrB1
% Mn* + ArrB2 -- kK5 w(n) --> Mn*.ArrB2

% case 0 --- calcium imaging data
% case 1 --- electrophysiology data
% case 2 --- store prameters values you want in data.mat
max_chan = 0;
maxchan_time = 0;

switch dataset    
    case 0        
        load('comparetoephys.mat')
    case 1
        load('incompleteset.mat')
    case 2
        load('data.mat')
    otherwise        
end

%% SPECIES: X = [
%% X(1)           G.GDP
%% X(2)           Ga.GTP
%% X(3)           Gbg
%% X(4)           PLC
%% X(5)           PLC*.Ga.GTP
%% X(6)           SecM      
%% X(7)           Channel-: closed channel
%% X(8)           SecM.Channel+: open channel
%% X(9)           PLC.Ga.GDP
%% X(10)          Ga.GDP    ];


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


t = 0;
tmax = 100;
counter =1; % counter counts the number of time iterations.
maxcounter=10000000;
tic;

K = [ kG1,      kG2,        kG3,    kG4*GTP,        kG5, ...    
       kP,      kI1,    kS*PIP2,         kO,         kC, ...
   kk1*Ki,      kk2,    kk3*ATP,  kk4*Arrb1,  kk5*Arrb2, ...
     kmax,       KM,        kI2,        kI3 ];
%%  K(1),    K(2),     K(3),      K(4),     K(5), 
%%  K(6),    K(7),     K(8),      K(9),    K(10),
%% K(11),   K(12),    K(13),     K(14),    K(15),  
%% K(16),   K(17),    K(18),     K(19)

no_rxns = 75;                   % number of reactions (total)
h = zeros(no_rxns,1);           % initialize the hazard vector
reaction_count = zeros(no_rxns,1); % it counts the number of each reactions' occuring.


%% store time, molecule numbers in every 'time_step' sec
time_step=0.25; % You record t, M, X in every time interval equal to time_step
tstore = zeros(floor(tmax/time_step)+1,1);
Xstore = zeros(floor(tmax/time_step)+1,size(X,2));
Mstore = zeros(floor(tmax/time_step)+1,size(M,2));
ttstore = zeros(floor(tmax/time_step)+1,1);

tstore(1,1) = t;
Xstore(1,:) = X;
Mstore(1,:) = M;
ttstore(1,1) = 0;
prev_t_index = 1; % it stores the previous time intex
%% store time, molecule numbers in every 'time_step' sec


%% Build function to account for increase in arrestin binding affinity
%% with more phosphates bound to melanopsin carboxyl tail.

W = @(n) 1-exp(-n*10); %it is alway 1

%% Build function to account for decrease in G-protein activation
%% with more phosphates bound to melanopsin carboxyl tail.

Y = @(n) exp(-n);

%% Begin the algorithm

for counter=1:maxcounter
    % check if the final time has been reached or exceeded
    if t>tmax
        break;
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
    h(22) = K(2)*M(28);
    h(27) = K(2)*M(35);
    h(32) = K(2)*M(42);

    % M0*.G.GDP -- kG3 --> M0*.G + GDP
    h(3)  = K(3)*M(2);  
    % And now for the other phosphorylated species
    h(8)  = K(3)*M(7);
    h(13) = K(3)*M(14);
    h(18) = K(3)*M(21);
    h(23) = K(3)*M(28);
    h(28) = K(3)*M(35);
    h(33) = K(3)*M(42);

    % M0*.G + GTP -- kG4 --> M0*.G.GTP
    h(4)  = K(4)*M(3);        
    % And now for the other phosphorylated species
    h(9)  = K(4)*M(8);
    h(14) = K(4)*M(15);
    h(19) = K(4)*M(22);
    h(24) = K(4)*M(29);
    h(29) = K(4)*M(36);
    h(34) = K(4)*M(43);

    % M0*.G.GTP -- kG5 --> M0* + Ga.GTP + Gbg
    h(5)  = K(5)*M(4);         
    % And now for the other phosphorylated species
    h(10) = K(5)*M(9);
    h(15) = K(5)*M(16);
    h(20) = K(5)*M(23);
    h(25) = K(5)*M(30);
    h(30) = K(5)*M(37);
    h(35) = K(5)*M(44);
    
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
    
    h(41) = K(11)*M(1);             % M0* + K -- kk1 --> M0*K 
    h(42) = K(12)*M(5);             % M0*K -- kk2 --> M0* + K
    h(43) = K(13)*M(5);             % M0*K + ATP -- kk3 --> M1*K + ADP
    
    h(44) = K(11)*M(6);             % M1* + K -- kk1 --> M1*K 
    h(45) = K(12)*M(10);            % M1*K -- kk2 --> M1* + K
    h(46) = K(13)*M(10);            % M1*K + ATP -- kk3 --> M2*K + ADP
    
    h(47) = K(11)*M(13);            % M2* + K -- kk1 --> M2*K 
    h(48) = K(12)*M(17);            % M2*K -- kk2 --> M2* + K
    h(49) = K(13)*M(17);            % M2*K + ATP -- kk3 --> M3*K + ADP
    
    h(50) = K(11)*M(20);            % M3* + K -- kk1 --> M3*K 
    h(51) = K(12)*M(24);            % M3*K -- kk2 --> M3* + K
    h(52) = K(13)*M(24);            % M3*K + ATP -- kk3 --> M4*K + ADP
    
    h(53) = K(11)*M(27);            % M4* + K -- kk1 --> M4*K 
    h(54) = K(12)*M(31);            % M4*K -- kk2 --> M4* + K
    h(55) = K(13)*M(31);            % M4*K + ATP -- kk3 --> M5*K + ADP
    
    h(56) = K(11)*M(34);            % M5* + K -- kk1 --> M5*K 
    h(57) = K(12)*M(38);            % M5*K -- kk2 --> M5* + K
    h(58) = K(13)*M(38);            % M5*K + ATP -- kk3 --> M6*K + ADP
            
    h(59) = K(11)*M(41);            % M6* + K -- kk1 --> M6*K 
    h(60) = K(12)*M(45);            % M6*K -- kk2 --> M6* + K
    
%%
% Arrestin Binding
    
    h(61) = W(1)*K(14)*M(6);        % M1* + ArrB1 -- kk4*w(n) --> M1*ArrB1
    h(62) = W(1)*K(15)*M(6);        % M1* + ArrB2 -- kk5*w(n) --> M1*ArrB2
    
    h(63) = W(2)*K(14)*M(13);       % M2* + ArrB1 -- kk4*w(n) --> M2*ArrB1
    h(64) = W(2)*K(15)*M(13);       % M2* + ArrB2 -- kk5*w(n) --> M2*ArrB2

    h(65) = W(3)*K(14)*M(20);       % M3* + ArrB1 -- kk4*w(n) --> M3*ArrB1
    h(66) = W(3)*K(15)*M(20);       % M3* + ArrB2 -- kk5*w(n) --> M3*ArrB2
    
    h(67) = W(4)*K(14)*M(27);       % M4* + ArrB1 -- kk4*w(n) --> M4*ArrB1
    h(68) = W(4)*K(15)*M(27);       % M4* + ArrB2 -- kk5*w(n) --> M4*ArrB2
    
    h(69) = W(5)*K(14)*M(34);       % M5* + ArrB1 -- kk4*w(n) --> M5*ArrB1
    h(70) = W(5)*K(15)*M(34);       % M5* + ArrB2 -- kk5*w(n) --> M5*ArrB2
    
    h(71) = W(6)*K(14)*M(41);       % M6* + ArrB1 -- kk4*w(n) --> M6*ArrB1
    h(72) = W(6)*K(15)*M(41);       % M6* + ArrB2 -- kk5*w(n) --> M6*ArrB2
        
%%
% SecM degradation
    h(73) = kmax*X(6)/(X(6)+KM);    % SecM -- delta --> 0, delta = kmax*SecM/(SecM+KM)

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
        reaction_count(1,1) = reaction_count(1,1) + 1;
    elseif hw(1) < r && r <= sum(hw(1:2))
        M(2) = M(2) - 1;            
        M(1) = M(1) + 1;           
        X(1) = X(1) + 1;    
        reaction_count(2,1) = reaction_count(2,1) + 1;
    elseif sum(hw(1:2)) < r && r <= sum(hw(1:3))
        M(2) = M(2) - 1;           
        M(3) = M(3) + 1;    
        reaction_count(3,1) = reaction_count(3,1) + 1;
    elseif sum(hw(1:3)) < r && r <= sum(hw(1:4))
        M(3) = M(3) - 1;           
        M(4) = M(4) + 1;  
        reaction_count(4,1) = reaction_count(4,1) + 1;
    elseif sum(hw(1:4)) < r && r <= sum(hw(1:5))
        M(4) = M(4) - 1;
        M(1) = M(1) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        reaction_count(5,1) = reaction_count(5,1) + 1;
        
        % using M1
    elseif sum(hw(1:5)) < r && r <= sum(hw(1:6))
        M(6) = M(6) - 1;
        X(1) = X(1) - 1;
        M(7) = M(7) + 1;
        reaction_count(6,1) = reaction_count(6,1) + 1;
    elseif sum(hw(1:6)) < r && r <= sum(hw(1:7))
        M(7) = M(7) - 1;
        M(6) = M(6) + 1;
        X(1) = X(1) + 1;
        reaction_count(7,1) = reaction_count(7,1) + 1;
    elseif sum(hw(1:7)) < r && r <= sum(hw(1:8))
        M(7) = M(7) - 1;
        M(8) = M(8) + 1;
        reaction_count(8,1) = reaction_count(8,1) + 1;
    elseif sum(hw(1:8)) < r && r <= sum(hw(1:9))
        M(8) = M(8) - 1;
        M(9) = M(9) + 1;
        reaction_count(9,1) = reaction_count(9,1) + 1;
    elseif sum(hw(1:9)) < r && r <= sum(hw(1:10))
        M(9) = M(9) - 1;
        M(6) = M(6) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        reaction_count(10,1) = reaction_count(10,1) + 1;
        
        % using M2
    elseif sum(hw(1:10)) < r && r <= sum(hw(1:11))
        M(13) = M(13) - 1;
        X(1) = X(1) - 1;
        M(14) = M(14) + 1;
        reaction_count(11,1) = reaction_count(11,1) + 1;
    elseif sum(hw(1:11)) < r && r <= sum(hw(1:12))
        M(14) = M(14) - 1;
        M(13) = M(13) + 1;
        X(1) = X(1) + 1;
        reaction_count(12,1) = reaction_count(12,1) + 1;
    elseif sum(hw(1:12)) < r && r <= sum(hw(1:13))
        M(14) = M(14) - 1;
        M(15) = M(15) + 1;
        reaction_count(13,1) = reaction_count(13,1) + 1;
    elseif sum(hw(1:13)) < r && r <= sum(hw(1:14))
        M(15) = M(15) - 1;
        M(16) = M(16) + 1;
        reaction_count(14,1) = reaction_count(14,1) + 1;
    elseif sum(hw(1:14)) < r && r <= sum(hw(1:15))
        M(16) = M(16) - 1;
        M(13) = M(13) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        reaction_count(15,1) = reaction_count(15,1) + 1;
        
        % using M3
    elseif sum(hw(1:15)) < r && r <= sum(hw(1:16))
        M(20) = M(20) - 1;
        X(1) = X(1) - 1;
        M(21) = M(21) + 1;
        reaction_count(16,1) = reaction_count(16,1) + 1;
    elseif sum(hw(1:16)) < r && r <= sum(hw(1:17))
        M(21) = M(21) - 1;
        M(20) = M(20) + 1;
        X(1) = X(1) + 1;
        reaction_count(17,1) = reaction_count(17,1) + 1;
    elseif sum(hw(1:17)) < r && r <= sum(hw(1:18))
        M(21) = M(21) - 1;
        M(22) = M(22) + 1;
        reaction_count(18,1) = reaction_count(18,1) + 1;
    elseif sum(hw(1:18)) < r && r <= sum(hw(1:19))
        M(22) = M(22) - 1;
        M(23) = M(23) + 1;
        reaction_count(19,1) = reaction_count(19,1) + 1;
    elseif sum(hw(1:19)) < r && r <= sum(hw(1:20))
        M(23) = M(23) - 1;
        M(20) = M(20) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        reaction_count(20,1) = reaction_count(20,1) + 1;
        
        % using M4
    elseif sum(hw(1:20)) < r && r <= sum(hw(1:21))
        M(27) = M(27) - 1;
        X(1) = X(1) - 1;
        M(28) = M(28) + 1;
        reaction_count(21,1) = reaction_count(21,1) + 1;
    elseif sum(hw(1:21)) < r && r <= sum(hw(1:22))
        M(28) = M(28) - 1;
        M(27) = M(27) + 1;
        X(1) = X(1) + 1;
        reaction_count(22,1) = reaction_count(22,1) + 1;
    elseif sum(hw(1:22)) < r && r <= sum(hw(1:23))
        M(28) = M(28) - 1;
        M(29) = M(29) + 1;
        reaction_count(23,1) = reaction_count(23,1) + 1;
    elseif sum(hw(1:23)) < r && r <= sum(hw(1:24))
        M(29) = M(29) - 1;
        M(30) = M(30) + 1;
        reaction_count(24,1) = reaction_count(24,1) + 1;
    elseif sum(hw(1:24)) < r && r <= sum(hw(1:25))
        M(30) = M(30) - 1;
        M(27) = M(27) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        reaction_count(25,1) = reaction_count(25,1) + 1;
        
        % using M5
    elseif sum(hw(1:25)) < r && r <= sum(hw(1:26))
        M(34) = M(34) - 1;
        X(1) = X(1) - 1;
        M(35) = M(35) + 1;
        reaction_count(26,1) = reaction_count(26,1) + 1;
    elseif sum(hw(1:26)) < r && r <= sum(hw(1:27))
        M(35) = M(35) - 1;
        M(34) = M(34) + 1;
        X(1) = X(1) + 1;
        reaction_count(27,1) = reaction_count(27,1) + 1;
    elseif sum(hw(1:27)) < r && r <= sum(hw(1:28))
        M(35) = M(35) - 1;
        M(36) = M(36) + 1;
        reaction_count(28,1) = reaction_count(28,1) + 1;
    elseif sum(hw(1:28)) < r && r <= sum(hw(1:29))
        M(36) = M(36) - 1;
        M(37) = M(37) + 1;
        reaction_count(29,1) = reaction_count(29,1) + 1;
    elseif sum(hw(1:29)) < r && r <= sum(hw(1:30))
        M(37) = M(37) - 1;
        M(34) = M(34) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        reaction_count(30,1) = reaction_count(30,1) + 1;
        
        % using M6
    elseif sum(hw(1:30)) < r && r <= sum(hw(1:31))
        M(41) = M(41) - 1;
        X(1) = X(1) - 1;
        M(42) = M(42) + 1;
        reaction_count(31,1) = reaction_count(31,1) + 1;
    elseif sum(hw(1:31)) < r && r <= sum(hw(1:32))
        M(42) = M(42) - 1;
        M(41) = M(41) + 1;
        X(1) = X(1) + 1;
        reaction_count(32,1) = reaction_count(32,1) + 1;
    elseif sum(hw(1:32)) < r && r <= sum(hw(1:33))
        M(42) = M(42) - 1;
        M(43) = M(43) + 1;
        reaction_count(33,1) = reaction_count(33,1) + 1;
    elseif sum(hw(1:33)) < r && r <= sum(hw(1:34))
        M(43) = M(43) - 1;
        M(44) = M(44) + 1;
        reaction_count(34,1) = reaction_count(34,1) + 1;
    elseif sum(hw(1:34)) < r && r <= sum(hw(1:35))
        M(44) = M(44) - 1;
        M(41) = M(41) + 1;
        X(2) = X(2) + 1;
        X(3) = X(3) + 1;
        reaction_count(35,1) = reaction_count(35,1) + 1;
        
        % PLC and G-protein activation/deactivation
        
    elseif sum(hw(1:35)) < r && r <= sum(hw(1:36))
        X(2) = X(2) - 1;
        X(4) = X(4) - 1;
        X(5) = X(5) + 1;
        reaction_count(36,1) = reaction_count(36,1) + 1;
    elseif sum(hw(1:36)) < r && r <= sum(hw(1:37))
        X(5) = X(5) - 1;
        X(9) = X(9) + 1;
        reaction_count(37,1) = reaction_count(37,1) + 1;
    elseif sum(hw(1:73)) < r && r <= sum(hw(1:74))
        X(9) = X(9) - 1;
        X(4) = X(4) + 1;
        X(10) = X(10) + 1;
        reaction_count(74,1) = reaction_count(74,1) + 1;
    elseif sum(hw(1:74)) < r && r <= sum(hw(1:75))
        X(10) = X(10) - 1;
        X(3) = X(3) - 1;
        X(1) = X(1) + 1;
        reaction_count(75,1) = reaction_count(75,1) + 1;
        
        % Second Messenger Creation
        
    elseif sum(hw(1:37)) < r && r <= sum(hw(1:38))
        X(6) = X(6) + 1;
        reaction_count(38,1) = reaction_count(38,1) + 1;
        
        % Channel Opening/Closing
        
    elseif sum(hw(1:38)) < r && r <= sum(hw(1:39))
        X(6) = X(6) - 1;
        X(7) = X(7) - 1;
        X(8) = X(8) + 1;
        reaction_count(39,1) = reaction_count(39,1) + 1;
    elseif sum(hw(1:39)) < r && r <= sum(hw(1:40))
        X(8) = X(8) - 1;
        X(6) = X(6) + 1;
        X(7) = X(7) + 1;
        reaction_count(40,1) = reaction_count(40,1) + 1;
        
        %Kinase Phosphorylation M0
        
    elseif sum(hw(1:40)) < r && r <= sum(hw(1:41))
        M(1) = M(1) - 1;
        M(5) = M(5) + 1;
        reaction_count(41,1) = reaction_count(41,1) + 1;
    elseif sum(hw(1:41)) < r && r <= sum(hw(1:42))
        M(5) = M(5) - 1;
        M(1) = M(1) + 1;
        reaction_count(42,1) = reaction_count(42,1) + 1;
    elseif sum(hw(1:42)) < r && r <= sum(hw(1:43))
        M(5) = M(5) - 1;
        M(10) = M(10) + 1;
        reaction_count(43,1) = reaction_count(43,1) + 1;
        
        % M1
        
    elseif sum(hw(1:43)) < r && r <= sum(hw(1:44))
        M(6) = M(6) - 1;
        M(10) = M(10) + 1;
        reaction_count(44,1) = reaction_count(44,1) + 1;
    elseif sum(hw(1:44)) < r && r <= sum(hw(1:45))
        M(10) = M(10) - 1;
        M(6) = M(6) + 1;
        reaction_count(45,1) = reaction_count(45,1) + 1;
    elseif sum(hw(1:45)) < r && r <= sum(hw(1:46))
        M(10) = M(10) - 1;
        M(17) = M(17) + 1;
        reaction_count(46,1) = reaction_count(46,1) + 1;
        
        %M2
        
    elseif sum(hw(1:46)) < r && r <= sum(hw(1:47))
        M(13) = M(13) - 1;
        M(17) = M(17) + 1;
        reaction_count(47,1) = reaction_count(47,1) + 1;
    elseif sum(hw(1:47)) < r && r <= sum(hw(1:48))
        M(17) = M(17) - 1;
        M(13) = M(13) + 1;
        reaction_count(48,1) = reaction_count(48,1) + 1;
    elseif sum(hw(1:48)) < r && r <= sum(hw(1:49))
        M(17) = M(17) - 1;
        M(24) = M(24) + 1;
        reaction_count(49,1) = reaction_count(49,1) + 1;
        
        %M3
        
    elseif sum(hw(1:49)) < r && r <= sum(hw(1:50))
        M(20) = M(20) - 1;
        M(24) = M(24) + 1;
        reaction_count(50,1) = reaction_count(50,1) + 1;
    elseif sum(hw(1:50)) < r && r <= sum(hw(1:51))
        M(24) = M(24) - 1;
        M(20) = M(20) + 1;
        reaction_count(51,1) = reaction_count(51,1) + 1;
    elseif sum(hw(1:51)) < r && r <= sum(hw(1:52))
        M(24) = M(24) - 1;
        M(31) = M(31) + 1;
        reaction_count(52,1) = reaction_count(52,1) + 1;
        
        %M4
        
    elseif sum(hw(1:52)) < r && r <= sum(hw(1:53))
        M(27) = M(27) - 1;
        M(31) = M(31) + 1;
        reaction_count(53,1) = reaction_count(53,1) + 1;
    elseif sum(hw(1:53)) < r && r <= sum(hw(1:54))
        M(31) = M(31) - 1;
        M(27) = M(27) + 1;
        reaction_count(54,1) = reaction_count(54,1) + 1;
    elseif sum(hw(1:54)) < r && r <= sum(hw(1:55))
        M(31) = M(31) - 1;
        M(38) = M(38) + 1;
        reaction_count(55,1) = reaction_count(55,1) + 1;
        
        %M5
        
    elseif sum(hw(1:55)) < r && r <= sum(hw(1:56))
        M(34) = M(34) - 1;
        M(38) = M(38) + 1;
        reaction_count(56,1) = reaction_count(56,1) + 1;
    elseif sum(hw(1:56)) < r && r <= sum(hw(1:57))
        M(38) = M(38) - 1;
        M(34) = M(34) + 1;
        reaction_count(57,1) = reaction_count(57,1) + 1;
    elseif sum(hw(1:57)) < r && r <= sum(hw(1:58))
        M(38) = M(38) - 1;
        M(45) = M(45) + 1;
        reaction_count(58,1) = reaction_count(58,1) + 1;
        
        %M6
        
    elseif sum(hw(1:58)) < r && r <= sum(hw(1:59))
        M(41) = M(41) - 1;
        M(45) = M(45) + 1;
        reaction_count(59,1) = reaction_count(59,1) + 1;
    elseif sum(hw(1:59)) < r && r <= sum(hw(1:60))
        M(45) = M(45) - 1;
        M(41) = M(41) + 1;
        reaction_count(60,1) = reaction_count(60,1) + 1;
        
        % Arrestin Binding
        % M1        
    elseif sum(hw(1:60)) < r && r <= sum(hw(1:61))
        M(6) = M(6) - 1;
        M(11) = M(11) + 1;
        reaction_count(61,1) = reaction_count(61,1) + 1;
    elseif sum(hw(1:61)) < r && r <= sum(hw(1:62))
        M(6) = M(6) - 1;
        M(12) = M(12) + 1;
        reaction_count(62,1) = reaction_count(62,1) + 1;
        
        %M2       
    elseif sum(hw(1:62)) < r && r <= sum(hw(1:63))
        M(13) = M(13) - 1;
        M(18) = M(18) + 1;
        reaction_count(63,1) = reaction_count(63,1) + 1;
    elseif sum(hw(1:63)) < r && r <= sum(hw(1:64))
        M(13) = M(13) - 1;
        M(19) = M(19) + 1;
        reaction_count(64,1) = reaction_count(64,1) + 1;
        
        %M3       
    elseif sum(hw(1:64)) < r && r <= sum(hw(1:65))
        M(20) = M(20) - 1;
        M(25) = M(25) + 1;
        reaction_count(65,1) = reaction_count(65,1) + 1;
    elseif sum(hw(1:65)) < r && r <= sum(hw(1:66))
        M(20) = M(20) - 1;
        M(26) = M(26) + 1;
        reaction_count(66,1) = reaction_count(66,1) + 1;
        
        %M4      
    elseif sum(hw(1:66)) < r && r <= sum(hw(1:67))
        M(27) = M(27) - 1;
        M(32) = M(32) + 1;
        reaction_count(67,1) = reaction_count(67,1) + 1;
    elseif sum(hw(1:67)) < r && r <= sum(hw(1:68))
        M(27) = M(27) - 1;
        M(33) = M(33) + 1;
        reaction_count(68,1) = reaction_count(68,1) + 1;
        
        %M5
    elseif sum(hw(1:68)) < r && r <= sum(hw(1:69))
        M(34) = M(34) - 1;
        M(39) = M(39) + 1;
        reaction_count(69,1) = reaction_count(69,1) + 1;
    elseif sum(hw(1:69)) < r && r <= sum(hw(1:70))
        M(34) = M(34) - 1;
        M(40) = M(40) + 1;
        reaction_count(70,1) = reaction_count(70,1) + 1;
        
        %M6
    elseif sum(hw(1:70)) < r && r <= sum(hw(1:71))
        M(41) = M(41) - 1;
        M(46) = M(46) + 1;
        reaction_count(71,1) = reaction_count(71,1) + 1;
    elseif sum(hw(1:71)) < r && r <= sum(hw(1:72))
        M(41) = M(41) - 1;
        M(47) = M(47) + 1;
        reaction_count(72,1) = reaction_count(72,1) + 1;
        
        % SecM degradation
    elseif sum(hw(1:72)) < r && r <= sum(hw(1:73))
        X(6) = X(6) - 1;
        reaction_count(73,1) = reaction_count(73,1) + 1;
    else
        M(1) = M(1);   % If the # of melanopsin cells is zero
    end
  
    
    %% store time, molecule numbers in every 'time_step' sec
    time_index = floor(t/time_step) + 1; % an index for the current time
    if time_index > floor(tmax/time_step)+1;
        time_index = floor(tmax/time_step)+1;
    end
    if time_index > prev_t_index
        for j = (prev_t_index+1):time_index
            tstore(j,1) = t;
          %% New
            if h_tot==0
                tstore(j,1) = (j-1)*time_step;
            end
          %% New
            Xstore(j,:) = X;
            Mstore(j,:) = M;
            ttstore(j,1) = tt;
        end
    end
    prev_t_index = time_index;    
    %% store time, molecule numbers in every 'time_step' sec

    
% %     if t>19.89
% %         keyboard;
% %     end
% %    pause(0.2)
% %     plot_matrx(tic,:) = [t X M];
% % figure(1)    
%     if the new t < tmax, repeat.  Else plot the results.
% %     subplot(2,1,1)
% %     plot(M)
% %     pause(0.1)
% %     subplot(2,1,2)
% %     plot(X)
% %    axis([0 14 0 100])
% %    pause(0.1)
% %    display(X(1))
% %    display(X(5))
% %    pause(0.1)

    %% record the maximal number of channels 
    if max_chan < X(8)
       max_chan = X(8);
       maxchan_time = t;
    end
    %% record the maximal number of channels
end
% % simstuff = {tstore,Mstore,Xstore};
% % save(sprintf('run%d',runnum),'simstuff')

save('results.mat','reaction_count','tstore','Mstore','Xstore','max_chan','maxchan_time')

%keyboard;
%% plotting
figure(1)
plot(tstore,Mstore(:,1)+Mstore(:,6)+Mstore(:,13)+Mstore(:,20)+Mstore(:,27)+Mstore(:,34)+Mstore(:,41));
%axis([0 tmax -1 100]);
xlabel('time (/sec)'); ylabel('# of cells');
figure(2)
plot(tstore, Xstore(:,8)./(Xstore(:,7)+Xstore(:,8))) 
% a ratio of the number of open channels out of the total number of
% channels
%% plotting 
toc