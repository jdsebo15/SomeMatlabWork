
clc;clear;clear globals;

run const.m
global Dpfr Lpfr RHOcat Init XSA P0;

Wmax= RHOcat*Lpfr*((Dpfr)^2)*(pi/4);                              %Determining the maximum catalyst weight [=] kg_cat
options = odeset('Reltol',1e-6,'Abstol',1e-6,'MaxStep',0.5);      %Options for the odesolver        

[W, Soln] = ode15s(@odes, [0 Wmax], Init, options);

W=W;                                                              %Weight of the Catalyst in the Reactor [=] kg_cat
Z=W./(XSA*RHOcat);                                                %Changing x axis to Length down the reactor [=] m

figure(1);
plot(Z,Soln(:,1:6));
title('Species Feed Rates vs. Length along the Reactor');
xlabel('Length (m)');
ylabel('Per Species Flow Rate (kmol/min)');
legend('A','B','C','D','E','I');
saveas(figure(1),'FeedRates.png');

figure(2);
plot(Z,Soln(:,8:9));
title('System Temperatures vs. Length along the Reactor');
xlabel('Length (m)');
ylabel('Temperature (K)');
legend('TR','THX');
saveas(figure(2),'Temps.png');

figure(3);
plot(Z,(Soln(:,7)*P0));
title('Pressure Drop vs. Length along the Reactor');
xlabel('Length (m)');
ylabel('Pressure Drop (Pa)');
legend('P');
saveas(figure(3),'Pressure.png');


