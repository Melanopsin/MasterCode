%% Melanopsin Stochastic Model
%% Sensitivity Study - Set up for LHS  1/06/14
%% Based on Model_LHS.m
%% results to be analyzed with PRCC and eFAST
clear all;
close all;

%% to set a random seed
seed_name = sum(100*clock);
rand('seed',seed_name)
%% to set a random seed

%% Sample size N
runs=1000;

%% LHS MATRIX  %%
Melanopsin_Parameter_settings_LHS;
%load('seed.mat')

    s_r = 0.5; % start ratio = 0.9
    e_r = 2; % end raio = 1.1

    GTP_varied = LHS_Call(s_r*GTP, GTP, e_r*GTP, 0 ,runs,'unif'); 
    PIP2_varied = LHS_Call(s_r*PIP2, PIP2, e_r*PIP2, 0 ,runs,'unif'); 
    Ki_varied = LHS_Call(s_r*Ki, Ki, e_r*Ki, 0 ,runs,'unif'); 
    ATP_varied = LHS_Call(s_r*ATP, ATP, e_r*ATP, 0 ,runs,'unif'); 
    kmax_varied = LHS_Call(s_r*kmax, kmax, e_r*kmax, 0 ,runs,'unif'); 
    KM_varied = LHS_Call(s_r*KM, KM, e_r*KM, 0 ,runs,'unif');
    
    kG1_varied = LHS_Call(s_r*kG1, kG1, e_r*kG1, 0 ,runs,'unif'); 
    kG2_varied = LHS_Call(s_r*kG2, kG2, e_r*kG2, 0 ,runs,'unif'); 
    kG3_varied = LHS_Call(s_r*kG3, kG3, e_r*kG3, 0 ,runs,'unif'); 
    kG4_varied = LHS_Call(s_r*kG4, kG4, e_r*kG4, 0 ,runs,'unif'); 
    kG5_varied = LHS_Call(s_r*kG5, kG5, e_r*kG5, 0 ,runs,'unif');
    
    kP_varied = LHS_Call(s_r*kP, kP, e_r*kP, 0 ,runs,'unif');
    kI1_varied = LHS_Call(s_r*kI1, kI1, e_r*kI1, 0 ,runs,'unif'); 
    kI2_varied = LHS_Call(s_r*kI2, kI2, e_r*kI2, 0 ,runs,'unif'); 
    kI3_varied = LHS_Call(s_r*kI3, kI3, e_r*kI3, 0 ,runs,'unif'); 
    kS_varied = LHS_Call(s_r*kS, kS, e_r*kS, 0 ,runs,'unif'); 
    kO_varied = LHS_Call(s_r*kO, kO, e_r*kO, 0 ,runs,'unif'); 
    kC_varied = LHS_Call(s_r*kC, kC, e_r*kC, 0 ,runs,'unif'); 
    
    kk1_varied = LHS_Call(s_r*kk1, kk1, e_r*kk1, 0 ,runs,'unif'); 
    kk2_varied = LHS_Call(s_r*kk2, kk2, e_r*kk2, 0 ,runs,'unif'); 
    kk3_varied = LHS_Call(s_r*kk3, kk3, e_r*kk3, 0 ,runs,'unif'); 
    kB1_varied = LHS_Call(s_r*kB1, kB1, e_r*kB1, 0 ,runs,'unif'); 
    kB2_varied = LHS_Call(s_r*kB2, kB2, e_r*kB2, 0 ,runs,'unif'); 
    kUB1_varied = LHS_Call(s_r*kUB1, kUB1, e_r*kUB1, 0 ,runs,'unif'); 
    kUB2_varied = LHS_Call(s_r*kUB2, kUB2, e_r*kUB2, 0 ,runs,'unif'); 
    kDe_varied = LHS_Call(s_r*kDe, kDe, e_r*kDe, 0 ,runs,'unif'); 
    

%% LHS MATRIX and PARAMETER LABELS
LHSmatrix=[GTP_varied PIP2_varied Ki_varied ATP_varied kmax_varied KM_varied...
           kG1_varied kG2_varied kG3_varied kG4_varied kG5_varied...
           kP_varied kI1_varied kI2_varied kI3_varied kS_varied kO_varied kC_varied...
           kk1_varied kk2_varied kk3_varied kB1_varied kB2_varied kUB1_varied kUB2_varied kDe_varied];   
       
for x=1:runs %Run solution 'runs' times choosing different values
    x
    [time,open_channels] = Melanopsin_Sensitivity(60,100,LHSmatrix,x);
    %% Melnopsin_Sensitivity.m (shot_mel.m) is a version of Melanopsin.m
    %% It gives one realization of stochastic simulation
    time_step = 0.25;
    open_channels_lhs(:,x)=open_channels(floor(time_points/time_step)+1);
    % under the assumption that open_channels is saved with
    % [time_points*species] matrix
end

%% Save the workspace
save Melanopsin_LHS.mat;

%% CALCULATE PRCC
alpha = 0.05;
[prcc sign sign_label]=PRCC(LHSmatrix,open_channels_lhs,1:length(time_points),PRCC_var,alpha);
PRCC_PLOT(LHSmatrix,open_channels_lhs,1:length(time_points),PRCC_var,y_var_label)