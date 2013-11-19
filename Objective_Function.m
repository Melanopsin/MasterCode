% The function sig_01 which has 3 inputs, pars, xx and yy and 2 outputs
% which are f and y0
function [f y0] = Objective_Function(pars, xx, yy)
% [f y0] = sig_ 01_f(pars, , y) Mean-squared errors from sigmoid curve. par
% Vector of parameters. xx,yy    Vectors of data.


%initcon refers to the initial conditions vector. The initial conditions of
%y1 through y15 are expressed in this matrix. par(1) is the y1 value,
%par(2) is the y12 value, and par(3) is the y15 value.
initco = [pars(1)^2 0 0 0 0 0 0 0 0.695486 0.113123 0 pars(2)^2 0 0 pars(3)^2];


%Given T and Y, the function returns the right side of the ODE as the
%vector yprime where xx are the times, and the next set in brackets are the
%initial conditions, intco ODE23s solves a stiff set of differential
%equations numerically
[T Y]=ode23s(@yprime, xx, initco);
y0 = Y(:,10);


% Calculates the error between yy and y0
f  = mean((yy - y0).^2);   % Means-squared errors.
    
function dy=yprime(~,y)


%rate constants for the differential equations listed below. The values for
%each rate constant is expressed in the par vector. The Parameters in order
%are:M*(y2),  M*.G.GDP(y3), M*.G(y4), M*.G.GTP(y5), Ga.GTP(y6),
%PLC*.Ga.GTP(y7), SecM(y8), Channel?(y9), SecM.Channel=(y10),
%PLC.Ga.GDP(y11), Ga.GDP(y13), Gby(y14),


KL=pars(4)^2;
KG1=pars(5)^2;
KG2=pars(6)^2;
KG3=pars(7)^2;
KG4=pars(8)^2;
KG5=pars(9)^2;
KP=pars(10)^2;
KS=pars(11)^2;
KC=pars(12)^2;
KI1=pars(13)^2;
KI2=pars(14)^2;
KI3=pars(15)^2;
%chemical constants. The values for g and i are in the par vector g = GTP
%and i=PIP2
g=pars(16)^2;
i=pars(17)^2;
%differential equations obtained from the chemical equations to represent
%the activation pathway of melanopsin.
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