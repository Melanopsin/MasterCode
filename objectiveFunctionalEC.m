
function output = objectiveFunctional(t,KnownData,IC,q)

%** ODE solver **
%[t,y] = ode45(@odeModel,t,IC,[ ],q);
[t,y] = GetAvg(t,IC,q);  %this is a function that has or gets the avg of the stochastic
%** Output **

ModelOutput = y(:,2); %choose the relevant column that we care abou & want to fit data
disp('check')
output = norm((KnownData - ModelOutput))^2;











