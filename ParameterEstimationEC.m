
function ParameterEstimationEC
%** Description of the program **
% This program is an implementation for estimating a set of parameters
% from the average of a stochastic output model
% via the least squares approach
%

CaImaging

clc       % Clear the display in the command window.
close all % Clear all stored variables from the Workspace.
clear all % Close all figures.

% Auxiliary parameters
MaxNumIter = 10; %Maximum number of iterations
TOL = 10^(-5); %Tolerance you want for your estimation

%** Model parameters **
load('calciumimaging.mat')

K = [kS kC kO kk3];

%** Independent variables **

tMin = 0.0;        % lower bound for the time
dt = 1;            % step size
tMax = 10;         % upper bound for the time; this should agree with the
t = [tMin:dt:tMax]';   % time vector t
n = length(t);         % Total number of points in t



%*********************** General Least Square *****************************

numIter = 0;  % counter

%K_opt = ones(length(K),1); %initial guess for optimum parameters, we can set to rand(p,1)
K_opt=K;

%** Auxiliary vectors **
convergence = [];

%% LOAD KNOWN DATA

KnownData_Exp = xlsread('brown_data.xls');
dt = 0.01; final_time = floor(15/dt); steps = 0.01;

KnownData = zeros(final_time,length(KnownData_Exp(1,:)));

for k = 1:final_time % last time point taken from experimental data set
    KnownData(k,:) = KnownData_Exp(max(find(KnownData_Exp(:,1)<steps)),:);
    steps = steps + dt;
end

% options=optimset('Algorithm','trust-region-reflective','tolfun',10^-6,'TolX',10^-6,'MaxFunEvals',1000,'MaxIter',1000);
% 'active-set', 'interior-point', 'sqp', 'trust-region-reflective'

% options=optimset('tolfun',10^-6,'TolX',10^-6,'MaxFunEvals',100,'MaxIter',100);

%while(numIter == 0 || ((norm(K - K_opt) >= TOL) && (numIter < MaxNumIter)))  % do / while

K = K_opt;  %K gets updated with the optimal current guess

disp('check 1')
[K_opt] = fminsearch(@Get_Error,KnownData,K);
disp('check 2')

kC = K_opt(1)
kS = K_opt(2);
kO = K_opt(3);
kk3 = K_opt(4);

save('calciumimaging.mat','kC','kS','kO','kk3','-append');

numIter = numIter + 1

convergence = [convergence; norm(K - K_opt)];

figure(1)
plot([1:1:numIter],convergence,'o');
xlabel('Number of iterations')
xlim([0 numIter]);
ylabel('Convergence')


%end % end while

%**************************************************************************

Values_Calc = Avg_Matrix();
t_Calc = Values_Calc(:,1);

%Y_Calc = sum(Values_Calc(:,2:8));

Y_Calc = Values_Calc(:,end);

figure(2)
plot(KnownData(:,1),KnownData(:,2),'b',t_Calc,Y_Calc,'m');

% disp('q_opt = '), disp(q_opt)
%disp('Parameter estimates:')
%disp(['beta = ',num2str(q(1))])
%disp(['gamma = ',num2str(q(2))])


