% outputs par2 to the workspace, so it can be manipulated after running the code
function [param_final] = Param_Estimation()

format long

% Initial Conditions

init_cond = [];

% Parameters

K = [];


%[T yCal] = ode23s(@yprime, [0 0.5], init_cond);

no_runs = 10;
no_species = 60;
final_time = 10; dt = 0.1; t_end = floor(final_time/dt); sizing = 1;

for j = 1:no_runs
    model_matrx = Melanopsin();
    for k = 1:final_time
        steps = dt;
        new_matrx(k,:) = model_matrx(max(find(model_matrx(:,1)<steps)),:);
        steps = steps + dt;
    end
    calc_matrx(sizing : t_end,:) = new_matrx;
    sizing = sizing + t_end;
end

avging_matrx = zeros(t_end,no_species);

for j = 1:no_runs
    avging_matrx_temp = calc_matrx((j-1)*t_end + j:j*t_end,:);
    avging_matrx = avging_matrx + avging_matrx_temp;
end

mean_matrx = avging_matrx/no_runs;
    

Calc_Results = [T yCal];


%Raw Electrophysiology Data

%Location of the data
filename='\\profile.ad.umbc.edu\Userfiles$\dthatch1\My Documents\data.xlsx';
sheet='Data Manipulation';
xrange='H1182:H1642';
yrange='L1182:L1642';

%Pulling the X and Y data values
times=xlsread(filename,sheet,xrange);
response=xlsread(filename,sheet,yrange);
plot(times,response,'.')
hold on


%fminunc
%adjustable options for the fminunc minimization functions
ops2 = optimset('LargeScale', 'off', 'Display', 'iter',  ...
    'TolX', 1e-8, 'TolFun', 1e-8);

Error  = mean((ydata - Calc_Results(:,end)).^2);

%outputs the final set of parameters  as a result of the minimization
param_final = fminunc([Error Calc_Results(:,end)], param_guess, ops2, times, response);


%plots the output of the minimization against the raw data
plot(times, response, '.', times, Calc_Results(:,end), 'r-');


%plot of model
% (:) means it looks at all the times
m = Calc_Results(:,end);
plot(T,m, '-k', 'LineWidth',2);
hold off


%plot specifications
legend('Data','model');
xlabel('t');
ylabel('solution');
axis([0 0.5 0 1])
axis on;
grid on;


    function dy=yprime(K,~,y)
        
        % DE's
        
        
        dy=[
            %1. dM/dt
            -K(1).*y(1);
            
            
            %2. dM*/dt
            -K(2)*y(2).*y(15)+K(3)*y(3)+K(4)*y(5)+K(1)*y(1);
            
            
            %3. dM*.G.GDP/dt
            K(2)*y(2).*y(15)-K(2)*y(3)-K(4)*y(3);
            
            
            %4. dM*.G/dt
            KG3.*y(3)-KG4.*y(4).*g;
            
            
            %5. dM*.G.GTP/dt
            KG4.*y(4).*g-KG5.*y(5);
            
            
            %6. dGa.GTP/dt
            KG5.*y(5)-KP.*y(12).*y(6);
            
            
            %7. dPLC*.Ga.GTP
            KP.*y(6).*y(12)-KI1.*y(7);
            
            
            %8. dSecM/dt
            K(8)*y(7)-KC.*y(8).*y(9);
            
            
            %9. dChannel?/dt
            -KC.*y(8).*y(9);
            
            
            %10. dSecM.Channel=/dt
            KC.*y(8).*y(9);
            
            
            %11. dPLC.Ga.GDP/dt
            KI1.*y(7)-KI2.*y(11);
            
            
            %12. dPLC/dt
            KI2.*y(11)-KP.*y(12).*y(6);
            
            
            %13. dGa.GDP/dt
            KI2.*y(11)-KI3.*y(13).*y(14);
            
            
            %14. dGby/dt
            KG5.*y(5)-KI3.*y(13).*y(14);
            
            
            %15. dG.GDP/dt
            KI3.*y(13).*y(14)+KG2.*y(3)-KG1.*y(2).*y(15)];
    end
end