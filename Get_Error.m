
function output = Get_Error(KnownData)

%** ODE solver **
%[t,y] = ode45(@odeModel,t,IC,[ ],q);
Model_Species = Avg_Matrix();  %this is a function that has or gets the avg of the stochastic
%** Output **

 %choose the relevant column that we care abou & want to fit data
 
ModelOutput = Model_Species(:,end);
disp('check')
 
output = norm((KnownData(:,2) - ModelOutput))^2;











