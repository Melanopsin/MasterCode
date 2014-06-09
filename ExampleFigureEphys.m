% % pathdef
set(0,'DefaultAxesLineStyleOrder','-|-.|--|:','DefaultLineLineWidth',3)
set(0,'DefaultTextFontSize',18)
set(0,'DefaultAxesFontSize',18)


fignum = 2
figure(fignum)
xlabel('Time (s)')
ylabel('Channel response (fraction of total channels)')
title('iprgc')
legend('Experiment','Model')
screen_size = get(0, 'ScreenSize');
set(fignum, 'Position', [0 0 0.75*screen_size(3) 0.5*screen_size(4) ] );
export_fig('EphysChannel','-pdf','-nocrop')

fignum = 3
figure(fignum)
title('Reaction count')
screen_size = get(0, 'ScreenSize');
set(fignum, 'Position', [0 0 0.75*screen_size(3) 0.5*screen_size(4) ] );
export_fig('EphysCount','-pdf','-nocrop')