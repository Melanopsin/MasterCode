%% Set initial conditions for Melanopsin.m, Melanopsin2.m,
%% and Melanopsin2_count.m
%% Initial values are stored in 'hwk.mat'

GTP = 1; 
PIP2 =1.5;
Ki = 1;
ATP = 1;
kmax =30; % increase -- increase deactivation, decrease -- decrease deactivation
KM =1;   % increase -- decrease deactivation, decrease -- increase deactivation


% kG1 = 6.13;  % increase -- nothing, decrease -- nothing 
% kG2 = .0106; % increase -- shift left, decrease -- shift right
% kG3 = .5759; % increase -- nothing clear, decrease -- decrease deactivation
% kG4 = .5443; % increase -- shfit left, decrease -- shift right, decrease deactivation
% kG5 = 3.5; %increase-- shift left & faster deactivation for t>40, decrease--shift right & slower deactivation for t>40 
% kP = 1;% increase --not clear *more runs, decrease --delays activation & changes channel dynamics
% kI1 = 15;% *increase --increase folds runs faster & faster deactivation equilibrium, decrease -- runs slower 
% kI2=.6730;% *increase -- shift slightly left, lowers deactivation curve but at same rate, decrease-- increase deactivation curve\ deactivation impossible 
% kI3=.00004;% increase -- shift right & slower activation , decrease -- shift to right, slow deactivation
% kS = 16.0733; %increase -- no clear difference, decrease -- shift left
% kO = 0.08; % increase -- shift right, decrease -- no clear difference 
% kC = 0.6; % increase -- slows simulation way down shift peak left , decrease -- shift peak up speeds up simulation
% kk1 = 50;% increase -- no clear change, decrease -- no clear change
% kk2 = .0861; % increase -- no clear change, decrease -- shift right
% kk3 = 10.0602; % increase -- no clear change, decrease -- no clear change
% kB1 = 9.7004; % increase -- shift right, decrease -- higher peak
% kB2 = 9.905; % increase -- no clear change, decrease -- no clear change
% 
% kUB1 = 10;      
% kUB2 = 10;       
% kDe = 10;

kG1 = 20;  % increase -- nothing, decrease -- nothing 
kG2 = 0.2; % increase -- shift left, decrease -- shift right
kG3 = 50.759; % increase -- nothing clear, decrease -- decrease deactivation
kG4 = 50.443; % increase -- shfit left, decrease -- shift right, decrease deactivation
kG5 = 50; %increase-- shift left & faster deactivation for t>40, decrease--shift right & slower deactivation for t>40 
kP = 10;% increase --not clear *more runs, decrease --delays activation & changes channel dynamics
kI1 = 150;% *increase --increase folds runs faster & faster deactivation equilibrium, decrease -- runs slower 
kI2=6.730;% *increase -- shift slightly left, lowers deactivation curve but at same rate, decrease-- increase deactivation curve\ deactivation impossible 
kI3=0.004;% increase -- shift right & slower activation , decrease -- shift to right, slow deactivation
kS = 50; %increase -- no clear difference, decrease -- shift left
kO = 0.79; % increase -- shift right, decrease -- no clear difference 
kC = 50.00; % increase -- slows simulation way down shift peak left , decrease -- shift peak up speeds up simulation
kk1 = 375.;% increase -- no clear change, decrease -- no clear change
kk2 = 80; % increase -- no clear change, decrease -- shift right
kk3 = 100.602; % increase -- no clear change, decrease -- no clear change
kB1 = 10%0.97004; % increase -- shift right, decrease -- higher peak
kB2 = 10%0.905; % increase -- no clear change, decrease -- no clear change

% % kG1 =  10;
% % kG2 =  0.2;
% % kG3 =  10;
% % kG4 =  7.6;
% % kG5 =  10.9;
% % kP =   5;
% % kI1 =  4.8;
% % kI2 =  7;
% % kI3 =  5;
% % kS  =  100;
% % kO  =  0.079;
% % kC  =  100;
% % 
% % %=======
% % kk1 =   15;
% % kk2 =   100;
% % kk3 =   100;
% % kB1 =   100;
% % kB2 =   200;

kUB1 = 0.01;      
kUB2 = 0.01;       
kDe = 20;


X = zeros(1,10);
M = zeros(1,47);

%% SPECIES: X = [
%% X(1)           G.GDP 
X(1) =floor(30); % critical for peak time
%% X(2)           Ga.GTP
%% X(3)           Gbg
%% X(4)           PLC
X(4) = floor(12);
%% X(5)           PLC*Ga.GTP
%% X(6)           SecM
%% X(7)           Channel-
X(7) =500;
%% X(8)           SecM.Channel+ ];
%% X(9)           PLC.Ga.GDP
%% X(10)          Ga.GDP
%% X(11)          Beta arrestin 1
X(11)=20;
%% X(12)          Beta arrestin 2
X(12)=20;


%% MELANOPSIN COMPLEXES: M = [
%% M(1)                        M0*
M(1) =floor(79);%100

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
M(11)=0;
%% M(12)                       M1*.ArrB2
M(12)=0;
%% M(13)                       M2*
%% M(14)                       M2*.G.GDP
%% M(15)                       M2*.G
%% M(16)                       M2*.G.GTP
%% M(17)                       M2*.K
%% M(18)                       M2*.ArrB1
M(18)=0;
%% M(19)                       M2*.ArrB2
%% M(20)                       M3*
%% M(21)                       M3*.G.GDP
%% M(22)                       M3*.G
%% M(23)                       M3*.G.GTP
%% M(24)                       M3*.K
%% M(25)                       M3*.ArrB1
M(25)=0;
%% M(26)                       M3*.ArrB2
%% M(27)                       M4*
%% M(28)                       M4*.G.GDP
%% M(29)                       M4*.G
%% M(30)                       M4*.G.GTP
%% M(31)                       M4*.K
%% M(32)                       M4*.ArrB1
M(32)=0;
%% M(33)                       M4*.ArrB2
%% M(34)                       M5*
%% M(35)                       M5*.G.GDP
%% M(36)                       M5*.G
%% M(37)                       M5*.G.GTP
%% M(38)                       M5*.K
%% M(39)                       M5*.ArrB1
M(39)=0;
%% M(40)                       M5*.ArrB2
%% M(41)                       M6*
%% M(42)                       M6*.G.GDP
%% M(43)                       M6*.G
%% M(44)                       M6*.G.GTP
%% M(45)                       M6*.K
%% M(46)                       M6*.ArrB1
M(46)=0;
%% M(47)                       M6*.ArrB2 ];
%% M(48)                       M0
%% M(49)                       MP
M(49)=0;

save('clhkah.mat');