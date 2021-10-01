function dSolnt = odes(W,Soln)

%M-file which compiles values for the derivatives
format long
global Dpfr RHOcat Void Dp U FHX Visc R dH298 ALPH MW P0 XSA;

Solnt = Soln';                                           %Transposed Matrix
FA=Soln(1);                                             %Feed Rate of A [=] kmol/min                                           
FB=Soln(2);                                             %Feed Rate of B [=] kmol/min
FC=Soln(3);                                             %Feed Rate of C [=] kmol/min
FD=Soln(4);                                             %Feed Rate of D [=] kmol/min
FE=Soln(5);                                             %Feed Rate of E [=] kmol/min
FI=Soln(6);                                             %Feed Rate of inerts [=] kmol/min
P=Soln(7)*P0;                                              %Inlet pressure [=] atm
TR=Soln(8);                                             %Inlet temperature of the reactor [=] K
THX=Soln(9);                                            %Inlet temperature of cooling fluid [=] K



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Mass Balance Equations to find species amounts                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=rateConstant(TR);                                     %Overall reaction rate constant; such that reaction rate [=] kmol/(kg*min)

i=1;
FNET=0;
while i<=(size(Solnt,2)-3)                               %Loop to sum molar flow rates [=] kmol/min
    F= Solnt(i);                                         %Total molar flow rate [=] kmol/min
    FNET = FNET + F;
    i=i+1;
end

    i=1;
    while i<=(size(Solnt,2)-3)                              %Loop to create a matrix of species compositions [=] (unitless)
        X(i)=Solnt(i)/FNET;
        i=i+1;
    end
        XA=X(1);                                                %Mole Fraction of A [=] (unitless)
        XB=X(2);                                                %Mole Fraction of B [=] (unitless)
        XC=X(3);                                                %Mole Fraction of C [=] (unitless)
        XD=X(4);                                                %Mole Fraction of D [=] (unitless)
        XE=X(5);                                                %Mole Fraction of E [=] (unitless)
        XI=X(6);                                                %Mole Fraction of I [=] (unitless)

CONT= (P/R/TR)/1000;                             %Concentration [=] kmol/m^3

i=1;
while i<=(size(Solnt,2)-3)
    SCON(i)= X(i)*CONT;           %Loop to create a matrix of concentrations [=] kmol/m^3
    i=i+1;
end
    CONA=SCON(1);                                           %Concentration of A [=] kmol/m^3
    CONB=SCON(2);                                           %Concentration of B [=] kmol/m^3
    CONC=SCON(3);                                           %Concentration of C [=] kmol/m^3
    COND=SCON(4);                                           %Concentration of D [=] kmol/m^3
    CONE=SCON(5);                                           %Concentration of E [=] kmol/m^3
    CONI=SCON(6);                                           %Concentration of I [=] kmol/m^3

i=1;
AvgMW=0;
while i<=(size(Solnt,2)-3)
    pMW=X(i)*MW(i);
    AvgMW= AvgMW + pMW;
    i=i+1;                                              %Average Molar weight [=] kg/kmol
end

%6 Material Balance Differential Equations; All have units [=]
%kmol/(kg_cat*min)
dSoln(1)= (k(2)*(CONC^ALPH(3)))-(k(1)*(CONA^ALPH(1))*(CONB^ALPH(2)))-(k(4)*(CONA^ALPH(5))*(COND^ALPH(6)));
dSoln(2)= (k(2)*(CONC^ALPH(3)))-(k(1)*(CONA^ALPH(1))*(CONB^ALPH(2)));
dSoln(3)= (k(1)*(CONA^ALPH(1))*(CONB^ALPH(2)))-(k(2)*(CONC^ALPH(3)))-(k(3)*(CONC^ALPH(4)));
dSoln(4)= (k(3)*(CONC^ALPH(4)))-(k(4)*(CONA^ALPH(5))*(COND^ALPH(6)));
dSoln(5)= (k(4)*(CONA^ALPH(5))*(COND^ALPH(6)));
dSoln(6)= 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Ergun Equation to Solve for Pressure Drop                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


AvgRHO= (P*AvgMW/R/TR)/1000;                     %Average Density [=] kg/m^3
XSA= (pi/4)*Dpfr^2;                                     %Cross sectional area of PFR [=] m^2
MFLOW= FNET*AvgMW;                                      %Mass flow rate [=] kg/min
G= MFLOW/XSA;                                           %Mass velocity [=] kg/(min*m^2)
Beta= (G/Dp)*((1-Void)/Void^3)*((150*(1-Void)*Visc/Dp)+1.75*G);   %((kg/(min*m^2))^2)/m

%Ergun Equation. Units are wonky, but they work out on paper [=] 1/(min^2*m)
dSoln(7)= -Beta/(AvgRHO*XSA*RHOcat*P0*3600);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Energy Balance to Solve for TR                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CP=heatCap(TR);                                         %Heat Capacity of each species in the reactor[=] kJ/(kmol*K)

CPW=heatCapW(THX);                                      %Heat Capacity of the coolant water [=] kJ/(kg*K)

ICP=intHeatCap(TR);                                     %Enthalpy departure from STP [=] kJ/kmol

    DeltaH(1)= dH298(1) + ICP(3) - ICP(2) - ICP(1);         %Change in Enthalpy for reaction 1 [=] kJ/kmol
    DeltaH(2)= dH298(2) + ICP(4) - ICP(3);                  %Change in Enthalpy for reaction 2 [=] kJ/kmol
    DeltaH(3)= dH298(3) + ICP(5) - ICP(4) - ICP(1);         %Change in Enthalpy for reaction 3 [=] kJ/kmol
    
    RRate(1)= (-k(1)*CONA^ALPH(1)*CONB^ALPH(2))+(k(2)*CONC^ALPH(3));    %Rate of reaction 1 [=] kmol/(kg_cat*min)
    RRate(2)= (-k(3)*CONC^ALPH(4));                                     %Rate of reaction 2 [=] kmol/(kg_cat*min)
    RRate(3)= (-k(4)*CONA^ALPH(5)*COND^ALPH(6));                        %Rate of reaction 3 [=] kmol/(kg_cat*min)
    
    ReacQR=0;
    i=1;
    while i<=size(dH298,2)
        dEnthSum=DeltaH(i)*RRate(i);
        ReacQR= ReacQR+dEnthSum;                          %Sum of Reaction Heats [=] kJ/(kg_cat*min)
        i=i+1;
    end
    
TempQR=4*U*(THX-TR)/(Dpfr*RHOcat)*(60/1000);             %Heat from temperature gradient [=] kJ/(kg_cat*min)

NetQR= TempQR-ReacQR;                                      %Sum of heat sources [=] kJ/(kg_cat*min)

i=1;
while i<=(size(Solnt,2)-3)
    MassFlow(i)=Solnt(i)*MW(i);                          %Mass Flow of each species [=] kg/min
    i=i+1;
end

NetCPR=0;
i=1;
while i<=(size(Solnt,2)-3)
    dNetCP=Solnt(i)*CP(i);
    NetCPR= NetCPR+dNetCP;                                %Sum of Feed*Cp [=] kJ/(min*K)
    i=i+1;
end

%Differential Equation describing temperature change of the PFR 
%[=] K/kg_cat
dSoln(8)=NetQR/NetCPR; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Energy Balance to Solve for THX                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

QHX= 4*U*(THX-TR)/(Dpfr*RHOcat)*(60/1000);               %Heat from temperature gradient [=] kJ/(kg_cat * min)
NetCPHX= FHX*CPW;                                        %Net heat capacity [=]  kJ/(min*K)

%Differential Equation describing the temperature change of the Heat
%Exchanger [=] K/kg_cat
dSoln(9)= -QHX/NetCPHX;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Making dSoln a column vector for ode15s
dSolnt=dSoln';


end

