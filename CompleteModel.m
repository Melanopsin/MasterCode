function CompleteModel
%{
 
 - light is now a "reactant" in the M0->M0* equation
 - can model series of different light intensities
 - see Abbey's guide on how to choose inputs
 
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
y(37) light
%}


% get conditions from user

t2 = input('\nDuration of first light level: ');
light_level1 = input('First light level: ');

length2 = input('\nDuration of second light level: ');
light_level2 = input('Second light level: ');

length3 = input('\nDuration of third light level: ');
light_level3 = input('Third light level: ');

interval4 = input('Add additional light interval? (1 for yes, 0 for no): ');
if interval4 == 1
    length4 = input('\nDuration of fourth light level: ');
    light_level4 = input('Fourth light level: ');
end

 
% dim backgroud light ------------------------
light_level = light_level1;
[T Y] = ode23s(@yprime, [0 t2], [1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 1 1 0 0 1 light_level]);
 
m = Y(:, 27);
plot(T,m,'-k','LineWidth',2);
hold on;


% bright flash ------------------------------
Ysize = size(Y); % find size of solution matrix
last_row = Ysize(1); % find index of last row (conditions at end of model)
new_initcons = Y(last_row, :); % set new initcons to conditions at end of last part
light_level = light_level2;
new_initcons(37) = light_level; % set initcon for light to new intensity
t1 = t2; % use end time of last part as starting time
t2 = t2 + length2; % set new ending time

[T Y] = ode23s(@yprime, [t1 t2], new_initcons);

m = Y(:, 27);
plot(T,m,'-r','LineWidth',2); % bright flash will be graphed in RED

% dim background again -----------------------
Ysize = size(Y);
last_row = Ysize(1);
new_initcons = Y(last_row, :);
light_level = light_level3;
new_initcons(37) = light_level;
t1 = t2;
t2 = t1 + length3;

[T Y] = ode23s(@yprime, [t1 t2], new_initcons);

m = Y(:, 27);
plot(T,m,'-k','LineWidth',2);



% fourth interval
if interval4 == 1
    Ysize = size(Y);
    last_row = Ysize(1);
    new_initcons = Y(last_row, :);
    light_level = light_level4;
    new_initcons(37) = light_level;
    t1 = t2;
    t2 = t1 + length4;

    [T Y] = ode23s(@yprime, [t1 t2], new_initcons);

    m = Y(:, 27);
    plot(T,m,'-r','LineWidth',2);
end
 
 
 
 
 
function dy=yprime(~,y)
 
 
% Constants
% g = GTP
g = 1;
% p = PIP2
p = 1;
 
kinase = 1;
arr1 = 1;
arr2 = 1;
 
% Michaelis Constants
kM = .02; % degradation by enzyme k
delta = 40; % max degradation rate
 
% Rate constants
kL = 5;
kG1 = 100;
kG2 = 0;
kG3 = 10;
kG4 = 10;
kG5 = 10;
kG6 = 10;
kG7 = 0;
kG8 = 1;
kG9 = 1;
kG10 = 1;
kG11 = 1;
kG12 = 0;
kG13 = 1;
kG14 = 1;
kG15 = 1;
kG16 = .01;
kG17 = 0;
kG18 = 1;
kG19 = 1;
kG20 = 1;
kP = 5;
kI1 = 12;
kI2 = 4;
kI3 = .2;
kS = 90;
kC1 = 100; % opening
kC2 = 6; % closing
kK1 = 10;
kK2 = 1;
kK3 = 1;
kK4 = 1;
kK5 = 1;
kK6 = 1;
kK7 = 1;
kK8 = 1;
kK9 = 1;
kK10 = 1;
kK11 = 1;
kK12 = 10;
kK13 = 10;
 
 
% Differential equations
dy=[
    % M0
    -kL.*y(1).*y(37)
    % M0*
    kL.*y(1).*y(37) - kG1.*y(2).*y(6) + kG2.*y(7) + kG5.*y(15) - kK1.*y(2).*kinase + kK2.*y(28);
    % M1*
    -kG6.*y(3).*y(6) + kG7.*y(8) + kG10.*y(16) + kK4.*y(29) - kK5.*y(3).*kinase;
    % M2*
    -kG11.*y(4).*y(6) + kG12.*y(9) + kG15.*y(17) + kK7.*y(30) - kK8.*y(4).*kinase;
    % M3*
    -kG16.*y(5).*y(6) + kG17.*y(10) + kG20.*y(18) + kK10.*y(31) - kK11.*y(5).*kinase - kK12.*y(5).*arr1 - kK13.*y(5).*arr2;
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
    kS.*p.*y(20) - kC1.*y(25).*y(26) - delta.*y(25)/(kM+y(25)).*y(25) + kC2.*y(27);
    % ClosedChannel
    -kC1.*y(25).*y(26) + kC2.*y(27);
    % SecM-OpenChannel
    kC1.*y(25).*y(26) - kC2.*y(27);
    % M0*-K
    kK1.*y(2).*kinase - kK2.*y(28) - kK3.*y(28);
    % M1*-K
    kK3.*y(28) - kK4.*y(29) + kK5.*y(3).*kinase - kK6.*y(29);
    % M2*-K
    kK6.*y(29) - kK7.*y(30) + kK8.*y(4).*kinase - kK9.*y(30);
    % M3*-K
    kK9.*y(30) - kK10.*y(31) + kK11.*y(5).*kinase;
    % ArrB1
    0
    % ArrB2
    0
    % M3-ArrB1
    kK12.*y(5).*arr1;
    % M3-ArrB2
    kK13.*y(5).*arr2;
    % K
    0
    % light
    0];
