
lanedata = textread('Lane.csv', '', 'delimiter', ',');

%open('Lane.csv')
plot(lanedata(5001:end,1)-1, -lanedata(5001:end,2)/max((abs(lanedata(5001:end,2)))))
set(gca, 'TickDir', 'out');

xlabel('time (s)', 'FontSize', 14,'FontWeight','Bold');

ylabel('Voltage', 'FontSize', 14,'FontWeight','Bold');

%h = legend('i = 1', 'i = 2');

%legend boxoff;

%set(h, 'FontSize', 14, 'Location', 'Best');

title('Electrophysiological recording from IPRGCs', 'FontSize', 18,'FontWeight','Bold');

axis([0 16 0 1]);

set(gca, 'TickDir', 'out');

set(gca,'XTick',0:4:16)

set(gca,'XTickLabel',{' 0', ' 4',' 8','12', '16'});

set(gca,'YTickLabel',{'  0','0.2','0.4','0.6','0.8','1.0'});

set(gca,'YTick',0:0.2:1)
