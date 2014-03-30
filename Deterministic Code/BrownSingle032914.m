function BrownSingle032914

%{

    This is the single flash response using the Brown electrophysiological
    data with final parameters chosen by fminunc.

    All melanopsin is instantaneously activated at time 0. The entire
    interval is in the dark.

    This model has been fit to data that has NOT been normalized. I
    subtracted the baseline (0.1) from the data and still need to refit the
    model.

y parameters
y(1)  M0
y(2)  M0*
y(3)  M1*
y(4)  M2*
y(5)  M3*
y(6)  G-GDP
y(7)  M0*-G-GDP
y(8)  M1*-G-GDP
y(9)  M2*-G-GDP
y(10) M3*-G-GDP
y(11) M0*-G
y(12) M1*-G
y(13) M2*-G
y(14) M3*-G
y(15) M0*-G-GTP
y(16) M1*-G-GTP
y(17) M2*-G-GTP
y(18) M3*-G-GTP
y(19) Galpha-GTP
y(20) PLC-Galpha-GTP
y(21) PLC
y(22) PLC-Galpha-GDP
y(23) Galpha-GDP
y(24) Gbetagamma
y(25) SecM
y(26) ClosedChannel
y(27) SecM-OpenChannel
y(28) M0*-K
y(29) M1*-K

y(30) M2*-K
y(31) M3*-K
y(32) ArrB1
y(33) ArrB2
y(34) M3-ArrB1
y(35) M3-ArrB2
y(36) K
y(37) light (not used)

%}


filename = '/afs/umbc.edu/users/a/b/abigail8/home/Documents/MATLAB/UBM/Data Plots/EphysFlash.xlsx';
sheet = 'Sheet1';
xrange = 'A1:A5001';
yrange = 'B1:B5001';

times = xlsread(filename, sheet, xrange);
response = xlsread(filename, sheet, yrange);
response = response - 0.1.*ones(size(response), 1);

% parameters in order (square rooted)
parsqrt = [
 9.946150532958260   0.096514576031749   2.899977459036370   2.856379277225519   2.837982301391146   3.159504477592329  -0.141074795249622   0.936924895047473   1.035022738089411   1.201658944930676 ...
 0.912900802391513  -0.078776637121018   1.122544725378410   1.011792784028624   1.011056169754815   0.185081303593623  -0.138572070612638   1.112410329378546   0.965748446866115   0.847701569133148 ...
 4.434910791386408   4.535447556550175   5.326093835636573  -0.604057412702856   8.709501673320611  10.043909792225392   2.140842606387945   3.034782690385838   1.187067268253176   1.255164969934004 ...
 1.206497412106647   1.034907328139204   1.114595558187959   1.095524386676172   1.036349384597015   0.942686354278032   1.093383358127027   1.079284462264678   3.142242488182647   0.140954071026512 ...
 0.403396570754774   3.133091091314507   1.000000000000000   1.000000000000000   1.000000000000000 .5 .01];

initcons = [0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 1 1 0 0 1 1 0];

[T Y] = ode23s(@yprime, times, initcons);
 
model = Y(:, 27);
plot(times, response, '.', times, model, 'r', 'Linewidth', 3);
xlabel('Time (seconds)');
ylabel('Normalized Response');
title('Single Flash Model Fit to Electrophysiological ipRGC Data');
whitebg('w');
axis on;
grid on;
legend('Data (normalized current)', 'Model (normalized open channel concentration)');


function dy=yprime(~,y)
 
 
% Constants
% g = GTP
g = 1;
% p = PIP2
p = 1;
 
kinase = 1;

% Michaelis Constants
kM = parsqrt(41)^2; % half-max SecM concentration for degradation by enzyme k
kmax = parsqrt(42)^2; % max degradation rate
 
% Rate constants
kL1 = parsqrt(45)^2;
kL2 = parsqrt(46)^2;
kG1 = parsqrt(1)^2;
kG2 = parsqrt(2)^2;
kG3 = parsqrt(3)^2;
kG4 = parsqrt(4)^2;
kG5 = parsqrt(5)^2;
kG6 = parsqrt(6)^2;
kG7 = parsqrt(7)^2;
kG8 = parsqrt(8)^2;
kG9 = parsqrt(9)^2;
kG10 = parsqrt(10)^2;
kG11 = parsqrt(11)^2;
kG12 = parsqrt(12)^2;
kG13 = parsqrt(13)^2;
kG14 = parsqrt(14)^2;
kG15 = parsqrt(15)^2;
kG16 = parsqrt(16)^2;
kG17 = parsqrt(17)^2;
kG18 = parsqrt(18)^2;
kG19 = parsqrt(19)^2;
kG20 = parsqrt(20)^2;
kP = parsqrt(21)^2;
kI1 = parsqrt(22)^2;
kI2 = parsqrt(23)^2;
kI3 = parsqrt(24)^2;
kS = parsqrt(25)^2;
kO = parsqrt(26)^2;
kC = parsqrt(27)^2;
kK1 = parsqrt(28)^2;
kK2 = parsqrt(29)^2;
kK3 = parsqrt(30)^2;
kK4 = parsqrt(31)^2;
kK5 = parsqrt(32)^2;
kK6 = parsqrt(33)^2;
kK7 = parsqrt(34)^2;
kK8 = parsqrt(35)^2;
kK9 = parsqrt(36)^2;
kK10 = parsqrt(37)^2;
kK11 = parsqrt(38)^2;
kB1 = parsqrt(39)^2;
kB2 = parsqrt(40)^2;
kUB1 = parsqrt(43)^2;
kUB2 = parsqrt(44)^2;
kDe = parsqrt(47)^2;
 
% Differential equations
dy=[
    % M0
    -kL1.*y(1).*y(37) + kDe.*y(38) + kL2.*y(38).*y(37);
    % M0*
    kL1.*y(1).*y(37) - kG1.*y(2).*y(6) + kG2.*y(7) + kG5.*y(15) - kK1.*y(2).*kinase + kK2.*y(28);
    % M1*
    -kG6.*y(3).*y(6) + kG7.*y(8) + kG10.*y(16) + kK4.*y(29) - kK5.*y(3).*kinase;
    % M2*
    -kG11.*y(4).*y(6) + kG12.*y(9) + kG15.*y(17) + kK7.*y(30) - kK8.*y(4).*kinase;
    % M3*
    -kG16.*y(5).*y(6) + kG17.*y(10) + kG20.*y(18) + kK10.*y(31) - kK11.*y(5).*kinase - kB1.*y(5).*y(32) - kB2.*y(5).*y(33);
    % G-GDP
    -kG1.*y(2).*y(6) + kG2.*y(7) - kG6.*y(3).*y(6) + kG7.*y(8) - kG11.*y(4).*y(6) + kG12.*y(9) - kG16.*y(5).*y(6) + kG17.*y(10) + kI3.*y(23).*y(24);
    % M0*-G-GDP
    kG1.*y(2).*y(6) - kG2.*y(7) - kG3.*y(7);
    % M1*-G-GDP
    kG6.*y(3).*y(6) - kG7.*y(8) - kG8.*y(8);
    % M2*-G-GDP
    kG11.*y(4).*y(6) - kG12.*y(9) - kG13.*y(9);
    % M3*-G-GDP
    kG16.*y(5).*y(6) - kG17.*y(10) - kG18.*y(10);
    % M0*-G
    kG3.*y(7) - kG4.*y(11).*g;
    % M1*-G
    kG8.*y(8) - kG9.*y(12).*g;
    % M2*-G
    kG13.*y(9) - kG14.*y(13).*g;
    % M3*-G
    kG18.*y(10) - kG19.*y(14).*g;
    % M0*-G-GTP
    kG4.*y(11).*g - kG5.*y(15);
    % M1*-G-GTP
    kG9.*y(12).*g - kG10.*y(16);
    % M2*-G-GTP
    kG14.*y(13).*g - kG15.*y(17);
    % M3*-G-GTP
    kG19.*y(14).*g - kG20.*y(18);
    % Galpha-GTP
    kG5.*y(15) - kP.*y(21).*y(19) + kG10.*y(16) + kG15.*y(17) + kG20.*y(18);
    % PLC*-Galpha-GTP
    kP.*y(21).*y(19) - kI1.*y(20);
    % PLC
    -kP.*y(21).*y(19) + kI2.*y(22);
    % PLC-Galpha-GDP
    kI1.*y(20) - kI2.*y(22);
    % Galpha-GDP
    kI2.*y(20) - kI3.*y(23).*y(24);
    % Gbetagamma
    kG5.*y(15) - kI3.*y(23).*y(24) + kG10.*y(16) + kG15.*y(17) + kG20.*y(18);
    % SecM
    kS.*p.*y(20) - kO.*y(25).*y(26) - kmax.*y(25)/(kM+y(25)) + kC.*y(27);
    % ClosedChannel
    -kO.*y(25).*y(26) + kC.*y(27);
    % SecM-OpenChannel
    kO.*y(25).*y(26) - kC.*y(27);
    % M0*-K
    kK1.*y(2).*kinase - kK2.*y(28) - kK3.*y(28);
    % M1*-K
    kK3.*y(28) - kK4.*y(29) + kK5.*y(3).*kinase - kK6.*y(29);
    % M2*-K
    kK6.*y(29) - kK7.*y(30) + kK8.*y(4).*kinase - kK9.*y(30);
    % M3*-K
    kK9.*y(30) - kK10.*y(31) + kK11.*y(5).*kinase;
    % ArrB1
    -kB1.*y(5).*y(32) + kUB1.*y(34)
    % ArrB2
    -kB2.*y(5).*y(33) + kUB2.*y(35)
    % M3-ArrB1
    kB1.*y(5).*y(32) - kUB1.*y(34);
    % M3-ArrB2
    kB2.*y(5).*y(33) - kUB2.*y(35);
    % K
    0;
    % light
    0;
    % M3
    kUB1.*y(34) + kUB2.*y(35) - kDe.*y(38) - kL2.*y(38).*y(37);
    ];

end
end







