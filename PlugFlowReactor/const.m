function [  ] = const(  )

%M file which holds all constant values for the system

global Dpfr Lpfr RHOcat Void Dp U FHX Visc dH298 ALPH MW Init P0;

Dpfr = .3;                                  %Diameter of the plug flow reactor, [=] m
Lpfr = 2;                                    %Length of the plug flow reactor, [=] m
RHOcat = 600;                                %Density of the catalyst [=] kg/m^3
Void = .65;                                  %Porosity of the catalyst bed [=] (unitless)
Dp = .005;                                   %Diameter of the catalyst particles [=] m
U = 1200;                                    %Net heat transfer coefficient [=] W/(m^2 * K)
FHX = 100;                                   %Outer jacket flow rate [=] kg/min
Visc = 0.00114;                              %Gas phase viscosity [=] kg/(m*min)
dH298 = [92770, -24164, 27289];              %Enthalpy at 298 [=] kJ/kmol
ALPH = [1.5, 1.2, 1.7, 1.2, 1.0, .9];        %Exponents for reaction rates [=] (unitless)
MW= [16, 120.12, 136.12, 136.12, 152.12, 28];%Molar Weight per species [=] kg/kmol
P0 = 40*101325;

%Initial conditions, Init(1-6)= Initial species? molar flow rate,
%Init(7)= Initial Pressure, Init(8)= Initial Temperature of the PFR bulk fluid
%Init(9)= Initial Temperature of the water in the heat exchanger

Init= [10;        %Inlet flow rate of A [=] kmol/min
        5;         %Inlet flow rate of B [=] kmol/min
        0;         %Inlet flow rate of C [=] kmol/min
        0;         %Inlet flow rate of D [=] kmol/min
        0;         %Inlet flow rate of E [=] kmol/min
        4;         %Inlet flow rate of inerts [=] kmol/min
        1;         %Inlet reduced pressure [=] dimensionless
        310;       %Inlet temperature of the reactor [=] K
        293;       %Inlet temperature of cooling fluid [=] K
        ];

end