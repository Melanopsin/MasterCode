% outputs par2 to the workspace, so it can be manipulated after running the code
function [par2]=model


%Given T and Y,  the
%parameters in order: M(y1), M*(y2),  M*.G.GDP(y3), M*.G(y4), M*.G.GTP(y5),
%Ga.GTP(y6), PLC*.Ga.GTP(y7), SecM(y8), Channel?(y9), SecM.Channel=(y10),
%PLC.Ga.GDP(y11), PLC(y12), Ga.GDP(y13), Gby(y14), G.GDP(y15)


%Given T and Y, the function returns the right side of the ODE as the
%vector yprime where
%[0 0.6] = 0 to 0.6 seconds, and the next set in brackets are the initial
%conditions
%ODE45 solves the functions numerically
%ODE23s solves a stiff set of differential equations numerically


%Shows more digits (10-15)
format long


%MODEL
%the list of variable optimization parameters, in order they are: M(y1), PLC(y12), G.GDP(y15) , %M*(y2),  M*.G.GDP(y3), M*.G(y4), M*.G.GTP(y5), %Ga.GTP(y6), PLC*.Ga.GTP(y7), SecM(y8), %Channel?(y9), SecM.Channel=(y10), %PLC.Ga.GDP(y11), Ga.GDP(y13), Gby(y14), 
par = [2.883 3.860 3.066 0.1 3.162 2.098 3.162 2.757 3.301 2.236 1830.332 0.089 2.191 2.646 2.236 3.847 3.564];


%takes the list of variable optimization parameters and squares them
paract = par.^2
%initcon refers to the initial conditions vector. The initial conditions of y1 through y15 are expressed in this matrix. par(1) is the y1 value, par(2) is the y12 value, and par(3) is the y15 value.
initco = [par(1)^2 0 0 0 0 0 0 0 0.695486 0.113123 0 par(2)^2 0 0 par(3)^2]
[T Y]=ode23s(@yprime, [0 0.5], initco);


%Raw Electrophysiology Data 
%Location of the data
filename='data.xlsx';
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


%outputs the final set of parameters  as a result of the minimization
par2 = fminunc(@sig_01_f, par, ops2, times, response)


% The function sig_01 which has 3 inputs, pars2, times and response and 2 outputs which are f2 and y02
[f2 y02] = sig_01_f(par2, times, response);


%plots the output of the minimization against the raw data
plot(times, response, '.', times, y02, 'r-');


%plot of model
% (:) means it looks at all the times
m = Y(:,10);
plot(T,m, '-k', 'LineWidth',2);
hold off


%plot specifications
legend('Data','model');
xlabel('t');
ylabel('solution');
axis([0 0.5 0 1])
axis on;
grid on;


function dy=yprime(~,y)


%rate constants for the differential equations listed below. The values for each rate constant is expressed in the par vector.
    KL=par(4)^2;
    KG1=par(5)^2;
    KG2=par(6)^2;
    KG3=par(7)^2;
    KG4=par(8)^2;
    KG5=par(9)^2;
    KP=par(10)^2;
    KS=par(11)^2;
    KC=par(12)^2;
    KI1=par(13)^2;
    KI2=par(14)^2;
    KI3=par(15)^2;


    %chemical constants. The values for g and i are in the par vector
    %g = GTP and i=PIP2
    g=par(16)^2;
    i=par(17)^2;


%differential equations obtained from the chemical equations to represent the activation pathway of melanopsin.
    dy=[
    %1. dM/dt
        -KL.*y(1);


    %2. dM*/dt
        -KG1.*y(2).*y(15)+KG2.*y(3)+KG5.*y(5)+KL.*y(1);


    %3. dM*.G.GDP/dt
        KG1.*y(2).*y(15)-KG2.*y(3)-KG3.*y(3);


    %4. dM*.G/dt
            KG3.*y(3)-KG4.*y(4).*g;


    %5. dM*.G.GTP/dt
        KG4.*y(4).*g-KG5.*y(5);


    %6. dGa.GTP/dt
        KG5.*y(5)-KP.*y(12).*y(6);


    %7. dPLC*.Ga.GTP
        KP.*y(6).*y(12)-KI1.*y(7);


    %8. dSecM/dt
        KS.*i.*y(7)-KC.*y(8).*y(9);


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