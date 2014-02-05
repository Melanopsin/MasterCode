clear all
%close all 

%clc

 set(0,'DefaultTextFontName','Arial')
 set(0,'DefaultAxesFontName','Arial')

no_runs = 1%500;
dataset = 0;

Melanopsin(no_runs,dataset); %Christy added this so you only have to run one script

no_species =9;
final_time =60; 
dt = 0.05; sizing = 1;


%new_matrx = zeros(size(model_matrx));

new_matrx = zeros(final_time,no_species);
calc_matrx = zeros(no_runs*final_time,no_species);

for j = 1:no_runs
    %model_matrx = Melanopsin2();
    %model_matrx = sort(10*rand(4));
    %disp('model_matrx = ')
    clear tout
    clear tout2
    clear Mout
    clear Mout2
    clear Xout
    clear Xout2
    clear model_matrx
    clear a
        a=load(sprintf('run%d.mat',j));
            tout = a.simstuff(:,1);
            tout2 = tout{1,1}';
            Mout = a.simstuff(:,2);
            Mout2 = Mout{1,1};
            Xout = a.simstuff(:,3);
            Xout2 = Xout{1,1};
            
            model_matrx = [tout2 Mout2(:,1) Mout2(:,6) Mout2(:,13) Mout2(:,20) Mout2(:,27) Mout2(:,34) Mout2(:,41) Xout2(:,8)];
            
    steps = dt;
    k=1;
    for counter = 1:dt:final_time
       
        new_matrx(k,:) = model_matrx(max(find(model_matrx(:,1)<steps)),:);
        steps = steps + dt;
        k=k+1;
    end
    
    calc_matrx(sizing:sizing - 1 + length(new_matrx(:,1)),:) = new_matrx;
    sizing = sizing + length(new_matrx(:,1));
end

t_end = length(new_matrx(:,1));

avging_matrx = zeros(t_end,no_species);

for j = 1:no_runs
    avging_matrx_temp = calc_matrx((j-1)*t_end + 1:j*t_end,:);
    avging_matrx = avging_matrx + avging_matrx_temp;
end

Avg_Matrix = avging_matrx/no_runs;
ActM=Avg_Matrix(:,2)+Avg_Matrix(:,3)+Avg_Matrix(:,4)+Avg_Matrix(:,5)+Avg_Matrix(:,6)+Avg_Matrix(:,7)+Avg_Matrix(:,8);
%figure(1)
%clf;
hold on
%%plot(0.5*Avg_Matrix(:,1),Avg_Matrix(:,9)/max(Avg_Matrix(:,9)),'g-','LineWidth',2);
Evansdata = textread('Evans1.csv', '', 'delimiter', ',');%Erika added this line
plot(Evansdata(1:60,1),Evansdata(1:60,2),'b-','LineWidth',2)
plot(Avg_Matrix(:,1),Avg_Matrix(:,9)/max(Avg_Matrix(:,9)),'r-','LineWidth',2);
legend('Calcium Imagining','Model');
title('kG1=57, kG2=kG3=kG4=10,kG5=150, kP=.009, kI1=200, kI2=kI3=10, kS=125, k0=.3, kC=2, kk1=5,kk2=2,kk3=kk4=kk5=10','FontSize', 12);
%set(gca, 'TickDir', 'out');

%xlabel('time (s)', 'FontSize', 14,'FontWeight','Bold');%Erika commented out 

%ylabel('Mean Number of Open Channels', 'FontSize', 14,'FontWeight','Bold');%Erika commented out 

%h = legend('Exp','Model');

%legend boxoff;

%set(h, 'FontSize', 14, 'Location', 'Best');

%title('Stochastic Model of Melanopsin Phototransduction Cascade', 'FontSize', 18,'FontWeight','Bold');

axis([0 60 0 1.1]);

%set(gca, 'TickDir', 'out');

%set(gca,'XTick',0:5:30)

%set(gca,'XTickLabel',{' 0', ' 5','10','15','20','25','30'});

%set(gca,'YTickLabel',{'  0','0.2','0.4','0.6','0.8','1'});

%set(gca,'YTick',0:0.2:1)
%hold off

xlabel('time (s)', 'FontSize', 14);

ylabel('Normalize # of Open Channels', 'FontSize', 14);

%h = legend('i = 9', 'i = 10');

%legend boxoff;

%set(h, 'FontSize', 14, 'Location', 'Best');
% 
% figure(2)
% hold on
% plot(Avg_Matrix(:,1),ActM/max(ActM),'r-','LineWidth',2);
% title('Total Melanopsin','FontSize', 14);
% %set(gca, 'TickDir', 'out');
    