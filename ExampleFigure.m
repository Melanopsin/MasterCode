pathdef
set(0,'DefaultAxesLineStyleOrder','-|-.|--|:','DefaultLineLineWidth',3)
set(0,'DefaultTextFontSize',18)
set(0,'DefaultAxesFontSize',18)

x=0:0.01:1;
y=sin(2*pi*x);
z=cos(2*pi*x);

fignum = 1
figure(fignum)
plot(x,y)
hold all
plot(x,z)
xlabel('I am the x-axis! (my units)')
ylabel('I am the y-axis! (my units)')
title('I am a supercool title')
legend('sin(x)','cos(x)')
screen_size = get(0, 'ScreenSize');
set(fignum, 'Position', [0 0 0.5*screen_size(3) 0.5*screen_size(4) ] );
export_fig('NoSoftMuFront','-pdf','-nocrop')