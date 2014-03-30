% This is a simulation of the data from Do and Yau 2013

function DoYauSim031014

fprintf('This model simulates the data from Do and Yau 2013:\nFlashes superimposed on varying levels of background light.\n\n');

% get number of trials (light levels)
num = input('How many different levels of background light? ');

% initialize a vector that will contain the light levels (L values)
dim_levels = zeros(num, 1);
% get light levels
for n = 1:num
    fprintf('Level %d: ', n);
    dim_levels(n) = input('');
end

% get the length of the dim light before the flash.
% This should be in the 1000 to 2000 range.
dim_length = input('How many seconds of dim light? ');

% get light level (L) and length of the bright superimposed flash
% level should be 100 - 1000
flash_level = input('Level of bright flash: ');
% length should be very small (especially for very bright flashes)
% length should be in range .001 to .5
flash_length = input('How long is the bright flash? (seconds) ');

% this vector contains the colors that the various trials will cycle
% through when plotting
color = ['r', 'y', 'g', 'c', 'b', 'm', 'k'];

% the variable we are plotting. 27 = open channel
plot_var = 27;

% initialize color of plot. This will correspond to the index of the color
% vector above.
col = 1;
for n = 1:num
    % run dim background
    % The last initcon is L
    [T Y] = ode23s(@yprime, [0 dim_length], [1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 1 1 0 0 1 dim_levels(n)]);
    
    Ysize = size(Y); % find size of solution matrix
    
    % get the final value of channel for normalization in next step
    % (this is the initial value for the bright flash, which we want to
    % normalize to 0)
    baseline = Y(Ysize(1), plot_var);
    
    last_row = Ysize(1); % index of last row of solution (conditions at end of model)
    new_initcons = Y(last_row, :); % set new initcons to conditions at end of last part
    new_initcons(37) = flash_level; % set initcon for light to flash level
    
    % run flash
    % first calculate end of this time interval
    t2 = dim_length + flash_length;
    [T Y] = ode23s(@yprime, [dim_length t2], new_initcons);
    Ysize = size(Y); % find size of solution matrix
    % get part of solution corresponding to channel
    y = Y(:, plot_var);
    % normalize by subtracting baseline
    normY = y - baseline.*ones(Ysize(1),1);
    
    % we want t=0 to correspond to the bright flash, so subtract the length
    % of the dim background
    times = T - dim_length.*ones(Ysize(1), 1);
    
    % I removed this plot to make the legend easier to set up. This part of
    % the response is negligible for short flash length
    
    %plot(times, normY, color(col), 'LineWidth', 2);
    
    % will use this later to set axis limits
    axismin = times(1);
    hold on;
    
    last_row = Ysize(1); % find index of last row (conditions at end of model)
    new_initcons = Y(last_row, :); % set new initcons to conditions at end of last part
    new_initcons(37) = dim_levels(n); % set initcon for light to new intensity
    
    % run more background
    % calculate new time interval
    t1 = t2;
    t2 = t2 + 10; % change 10 to larger number to see more of the response
    [T Y] = ode23s(@yprime, [t1 t2], new_initcons);
    Ysize = size(Y);
    y = Y(:, plot_var);
    
    % normalize the response y subtracting baseline
    normY = y - baseline.*ones(Ysize(1),1);
    
    Tsize = size(T);
    
    % calculate times to plot and maximum value on x axis
    times = T - dim_length.*ones(Tsize(1), 1);
    axismax = times(Tsize(1));
    
    % plot the dim background after the flash
    plot(times, normY, color(col), 'LineWidth', 2);
    
    % set axis
    axis([axismin, axismax, 0, 0.8]);
    hold on;
    
    % cycle through colors in color vector above
    col = col + 1;
    if col == 8
        col = 1;
    end
    
   
    
end    
    
    
    
function dy=yprime(~,y)
 
% parameters in order
% see parameter definitions below to see which is which
par = [
          98.9259104242659
       0.00931506338658826
          8.40986926291904
          8.15890257536338
          8.05414354300938
          9.98246854392598
        0.0199020978547228
         0.877828258959718
           1.0712720683621
          1.44398421993191
         0.833387875007068
       0.00620575855609655
          1.26010666047489
          1.02372463781239
          1.02223457839928
        0.0342550889399148
        0.0192022187538739
          1.23745674090808
         0.932670062624313
         0.718597950310801
          19.6684337275556
           20.570284538217
          28.3672755460059
         0.364885357841268
          75.8554193975745
          100.880123914361
          4.58320706532593
          9.20990597786551
          1.40912869935806
          1.57543910174943
          1.45563600542004
          1.07103317783623
          1.24232325833233
           1.2001736818022
          1.07402004695461
         0.888657562542007
          1.19548716782913
          1.16485495048596
          9.87368785454027
         0.019868050138947
         0.162728793296711
          9.81625978647433
                         1
                         1
                         1];

% Constants
% g = GTP
g = 1;
% p = PIP2
p = 1; 
kinase = 1;
arr1 = 1;
arr2 = 1;
 
% Michaelis Constants
kM = par(41); % half-maxrate concentration for degradation by enzyme
delta = par(42); % max degradation rate
 
% Rate constants
kL = par(45);
kG1 = par(1);
kG2 = par(2);
kG3 = par(3);
kG4 = par(4);
kG5 = par(5);
kG6 = par(6);
kG7 = par(7);
kG8 = par(8);
kG9 = par(9);
kG10 = par(10);
kG11 = par(11);
kG12 = par(12);
kG13 = par(13);
kG14 = par(14);
kG15 = par(15);
kG16 = par(16);
kG17 = par(17);
kG18 = par(18);
kG19 = par(19);
kG20 = par(12);
kP = par(21);
kI1 = par(22);
kI2 = par(23);
kI3 = par(24);
kS = par(25);
kC1 = par(26);
kC2 = par(27);
kK1 = par(28);
kK2 = par(29);
kK3 = par(30);
kK4 = par(31);
kK5 = par(32);
kK6 = par(33);
kK7 = par(34);
kK8 = par(35);
kK9 = par(36);
kK10 = par(37);
kK11 = par(38);
kK12 = par(39);
kK13 = par(40);
kR1 = par(43);
kR2 = par(44);
 
 
% Differential equations
dy=[
    % M0
    -kL.*y(1).*y(37) + kR1.*y(34) + kR2.*y(35);
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
    kS.*p.*y(20) - kC1.*y(25).*y(26) - delta.*y(25)/(kM+y(25)) + kC2.*y(27);
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
    kK12.*y(5).*arr1 - kR1.*y(34);
    % M3-ArrB2
    kK13.*y(5).*arr2 - kR2.*y(35);
    % K
    0
    % light
    0];
        