%% PARAMETER BASELINE VALUES and
%% INITIAL CONDITION FOR THE STOCHASTIC MODEL (y0)
data_set=6;
switch data_set    
    case 1
        load('ephys.mat')     % electrophysiology
    case 2
        load('multidata.mat') % multi flash
    case 3
        load('beta1.mat')     % overexpressed beta arrestin1
    case 4
        load('beta2.mat')     % overexpressed beta arrestin2        
    otherwise 
        Set_IC_caimage
        load('caimage.mat')   % calcium imaging is the default  
end


%% Parameter Labels (26 parameters)
PRCC_var={'GTP', 'PIP2', 'Ki', 'ATP', 'kmax', 'KM',...
           'kG1', 'kG2', 'kG3', 'kG4', 'kG5',...
           'kP', 'kI1', 'kI2', 'kI3', 'kS', 'kO', 'kC',...
           'kk1', 'kk2', 'kk3', 'kB1', 'kB2', 'kUB1', 'kUB2', 'kDe'};

%% Variables Labels
y_var_label={'Open Channels'};

%% TIME SPAN OF THE SIMULATION
t_end=60; % length of the simulations
%% Multiple time points are not working
% time_points=[1 10 20 30 40 50 t_end]; % time points of interest for the US analysis
%% this does not work
time_points=[t_end];