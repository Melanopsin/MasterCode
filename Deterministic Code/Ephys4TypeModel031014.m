function Ephys4TypeModel031014

% get user input for model type
type = input('\nWhat would you like to model?\nEnter 0 for single flash\nEnter 1 for flash superimposed on background light\nEnter 2 for multiple flashes in the dark\nEnter 3 for custom\n');

% variable that we are plotting
plotvar = 27;

% if single flash
if type == 0
    flashLevel = input('\nLight level: \n');
    flashLength = input('\nHow long is the light flash? (seconds)\n');
    darkLength = input('\nHow long is the dark interval after the flash? (seconds)\n');

    % flash
    [T Y] = ode23s(@yprime, [0 flashLength], [1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 1 1 0 0 1 flashLevel]);
    m = Y(:, 27);
    plot(T,m,'-r','LineWidth',2);
    hold on;
    
    % dark
    Ysize = size(Y); % find size of solution matrix
    last_row = Ysize(1); % find index of last row (conditions at end of model)
    new_initcons = Y(last_row, :); % set new initcons to conditions at end of last part
    new_initcons(37) = 0; % set initcon for light to new intensity
    [T Y] = ode23s(@yprime, [flashLength darkLength], new_initcons);
    m = Y(:, 27);
    plot(T,m,'-k','LineWidth',2);
    hold off;
    
elseif type == 1
    dimLevel = input('\nLight level of the dim background: \n');
    dimLength1 = input('\nHow long is the first dim background interval? (seconds)\n');
    flashLevel = input('\nLight level of the bright flash:\n');
    flashLength = input('\nHow long is the bright flash?\n');
    dimLength2 = input('\nHow long is the second dim background interval?\n');
    
    % dim
    [T Y] = ode23s(@yprime, [0 dimLength1], [1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 1 1 0 0 1 dimLevel]);
    m = Y(:, plotvar);
    plot(T,m,'-k','LineWidth',2);
    hold on;
    
    % flash
    Ysize = size(Y); % find size of solution matrix
    last_row = Ysize(1); % find index of last row (conditions at end of model)
    new_initcons = Y(last_row, :); % set new initcons to conditions at end of last part
    new_initcons(37) = flashLevel; % set initcon for light to new intensity
    t2 = dimLength1 + flashLength;
    [T Y] = ode23s(@yprime, [dimLength1 t2], new_initcons);
    m = Y(:, plotvar);
    plot(T,m,'-r','LineWidth',2);
    hold on;
    
    % dim
    Ysize = size(Y); % find size of solution matrix
    last_row = Ysize(1); % find index of last row (conditions at end of model)
    new_initcons = Y(last_row, :); % set new initcons to conditions at end of last part
    new_initcons(37) = dimLevel; % set initcon for light to new intensity
    t1 = t2;
    t2 = t1 + dimLength2;
    [T Y] = ode23s(@yprime, [t1 t2], new_initcons);
    m = Y(:, plotvar);
    plot(T,m,'-k','LineWidth',2);
    hold off;
    
elseif type == 2
    numFlashes = input('\nHow many flashes?\n');
    flashLength = input('\nHow long is each flash?\n');
    flashLevel = input('\nLight level of flashes:\n');
    darkLength = input('\nHow many seconds between flashes?\n');
    
    t1 = 0;
    t2 = flashLength;
    for n = 1:numFlashes
        % flash
        if n ~= 1
            Ysize = size(Y); % find size of solution matrix
            last_row = Ysize(1); % find index of last row (conditions at end of model)
            new_initcons = Y(last_row, :); % set new initcons to conditions at end of last part
            new_initcons(37) = flashLevel; % set initcon for light to new intensity
            [T Y] = ode23s(@yprime, [t1 t2], new_initcons);
        
        else
            [T Y] = ode23s(@yprime, [t1 t2], [1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 1 1 0 0 1 flashLevel]);
        end
        m = Y(:, plotvar);
        plot(T,m,'-r','LineWidth',2);
        hold on;
    
        t1 = t2;
        t2 = t1 + darkLength;
        % dark
        Ysize = size(Y); % find size of solution matrix
        last_row = Ysize(1); % find index of last row (conditions at end of model)
        new_initcons = Y(last_row, :); % set new initcons to conditions at end of last part
        new_initcons(37) = 0; % set initcon for light to new intensity
        [T Y] = ode23s(@yprime, [t1 t2], new_initcons);
     
        m = Y(:, plotvar);
        plot(T,m,'-k','LineWidth',2);
        hold on;
        
        t1 = t2;
        t2 = t1 + flashLength;
    end
    
   
elseif type == 3
    
    numIntervals = input('\nHow many intervals?\n');
    
    % matrix with flash level (column 1), length (column 2) for each
    % interval (row)
    conditions = zeros(numIntervals, 2);    
    for n = 1:numIntervals
        fprintf('\nFlash level %d:\n', n);
        conditions(n, 1) = input('');
        fprintf('\nFlash length %d:\n', n);
        conditions(n, 2) = input('');
    end
    
   
    % set time parameters for interval 1
    t1 = 0;
    t2 = conditions(n,2);
    for n = 1:numIntervals
        % reset initcons if not first interval
        if n ~= 1
            Ysize = size(Y); % find size of solution matrix
            last_row = Ysize(1); % find index of last row (conditions at end of model)
            new_initcons = Y(last_row, :); % set new initcons to conditions at end of last part
            new_initcons(37) = conditions(n,1); % set initcon for light to new intensity
            [T Y] = ode23s(@yprime, [t1 t2], new_initcons);
        
        else % default initcons with flash level 1
            [T Y] = ode23s(@yprime, [t1 t2], [1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 1 1 0 0 1 conditions(n, 1)]);
        end
        m = Y(:, plotvar);
        
        % color of plot will alternate between red and black
        if mod(n, 2) == 0
            color = 'k';
        else color = 'r';
        end
        
        
        plot(T, m, color,'LineWidth',2);
        hold on;
        
    
        % set new time parameters if not on last interval
        if n ~= numIntervals
            t1 = t2;
            t2 = t1 + conditions(n+1, 2);
        end
        
      end

else fprintf('Invalid input\n');
    
end
    
    
    
function dy=yprime(~,y)
 
 
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
kM = par(41); % degradation by enzyme k
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
        
        
        
        