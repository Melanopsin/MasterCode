%% Rates from electrophysiology experiments
%% This is an incomplete set of rate constants
Moley = 10;
arrb1 = 10*Moley;
arrb2 = 10*Moley;

kG1 = 10.0;
kG2 = 4.4;
kG3 = 10.0;
kG4 = 7.6*14.8;
kG5 = 10.9;
kp = 5.0;
kI = 4.8;
%kI2 = 7.0;
%kI3 = 5.0;
kS = 3796000*12.7;
kO = 0.0079;
kC = 0;
kk1 = 0000.000;
kk2 = 0000.000;
kk3 = 0000.000;
kk4 = 0000.000;
kk5 = 0000.000;

X = zeros(13,1);
M = zeros(47,1);

%% SPECIES: X = [
%% X(1)           G.GDP 
X(1) = 9.4*Moley;
%% X(2)           G.GTP 
%% X(3)           Ga.GTP
%% X(4)           Gbg
%% X(5)           PLC
X(5) = 14.9 * Moley;
%% X(6)           PLC*.Ga.GTP
%% X(7)           PLC.Ga.GDP
%% X(8)           Ga.GDP
%% X(9)           G.GDP
%% X(10)          PIP2
%% X(11)          SecM
%% X(12)          Channel-
X(12) = floor(0.703  * Moley);
%% X(13)          SecM.Channel+ ];
X(13) = floor(0.114 * Moley);

%% MELANOPSIN COMPLEXES: M = [
%% M(1)                        M0*
M(1) = 8.9*0.92*Moley;
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



save('incompleteset')
