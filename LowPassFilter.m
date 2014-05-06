Exp_data =dlmread('EphysGraph.csv'); % single flash 
Time = Exp_data(:,2);
Data_store(:,1) = Exp_data(:,1);
Data_store(:,2) = Exp_data(:,3);
Data_store(:,3) = Exp_data(:,4);
Data_store(:,4) = Exp_data(:,5);

%% low pass filter (filter frequencies above cut-off frequency)
Fs = 15001; % length of data
fc = 1; % cut-off frequency
Wn = (2/Fs)*fc;
b = fir1(20,Wn,'low',kaiser(21,3));
fvtool(b,1,'Fs',Fs);
y = filter(b,1,Data_store(:,1));
newdata(:,1)=Time;
newdata(:,2)=(y-0.1e-10)/(max(y)-0.1e-10);

%plot(Time,Data_store(:,1),'b');
% % hold on;
% % plot(Time,(y-0.1e-10)/(max(y)-0.1e-10),'r','linewidth',2);
% % legend('Original Signal','Filtered Data','Location','NorthEast');
% % xlabel('Seconds'); ylabel('Amplitude');
%% low pass filter