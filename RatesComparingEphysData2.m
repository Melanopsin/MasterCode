clear all
close all
clc

%% Data from calcium imaging experiments

%% These are compared to the electrophysiology data

Moley = 10;
arrb1 = floor(10*Moley);
arrb2 = floor(10*Moley);

kG1 = 10*5.7;
kG2 = 0;
kG3 = 12.1;
kG4 = 12.19;
kG5 = 12.3;
kp = 4.7;
%kI1 = 1;
%kI2 = 22;
kI = 1000*0.1;
kS = 3104.9;
kO = 30; % try to change from 30
kC = 0.3; % try to change from 10
kk1 = 0.25; % try to change from 0.25
kk2 = 16.0; % try to change from 16.0
kk3 = 10*14.0*0.3; % try to change from 10*14.0*0.3
kk4 = 0.0105; % try to change from 0.0105
kk5 = 0.0088; % try to change from 0.0088

X = zeros(1,13);
M = zeros(1,47);

%% SPECIES: X = [
%% X(1)           G.GDP 
X(1) = floor(3*Moley);
%% X(2)           G.GTP 
%% X(3)           Ga.GTP
%% X(4)           Gbg
%% X(5)           PLC
X(5) = floor(1.2*Moley);
%% X(6)           PLC*.Ga.GTP
%% X(7)           PLC.Ga.GDP
%% X(8)           Ga.GDP
%% X(9)           G.GDP
%% X(10)          PIP2
X(10) = floor(1);
%% X(11)          SecM
%% X(12)          Channel-
X(12) = floor(50*Moley);
%% X(13)          SecM.Channel+ ];

%% MELANOPSIN COMPLEXES: M = [
%% M(1)                        M0*
M(1) = floor(8.6*0.92*Moley);

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


save('comparetoephys2')
