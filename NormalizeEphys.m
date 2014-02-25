% % figure(2)
% % Edata = dlmread('EphysGraph.csv')
% % 
% % plot(Edata(1291:end,2)-1.29,Edata(1291:end,5),'r')
% % hold all
% % opchan=Mx(:,8)./(Mx(:,7)+Mx(:,8));
% % plot(tstore(:,1,1),opchan/max(opchan)) 

figure(3)
plot(opchan/max(opchan)) 
