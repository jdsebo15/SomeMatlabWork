function ICP = intHeatCap( TR )

%m-file which calculates the integrals of the per-species heat capacities
%in order to calculate DeltaH @ TR instead of 298K

    CPA = @(T) (5.5 + .228.*T - (1.57.*10^(-4)).*T.^2);         %Heat Capacity of A [=] kJ/(kmol*K)
    CPB = @(T) (23.03 + .18.*T - (1.04.*10^(-4)).*T.^2);        %Heat Capacity of B [=] kJ/(kmol*K)
    CPC = @(T) (4.04 + .21.*T - (1.21.*10^(-4)).*T.^2);         %Heat Capacity of C [=] kJ/(kmol*K)
    CPD = @(T) (4.40 + .192.*T - (1.35.*10^(-4)).*T.^2);        %Heat Capacity of D [=] kJ/(kmol*K)
    CPE = @(T) (7.32 + .20.*T - (1.25.*10^(-4)).*T.^2);         %Heat Capacity of E [=] kJ/(kmol*K)
    CPI = @(T) (30.25 - (8.79.*10^(-3).*T) + (1.87.*10^(-5).*T.^2) - (7.6.*10^(-9).*T.^3));        %Heat Capacity of I [=] kJ/(kmol*K)
    
    ICP(1)= integral(CPA,298,TR);                               %Enthalpy departure of A [=] kJ/kmol
    ICP(2)= integral(CPB,298,TR);                               %Enthalpy departure of B [=] kJ/kmol
    ICP(3)= integral(CPC,298,TR);                               %Enthalpy departure of C [=] kJ/kmol
    ICP(4)= integral(CPD,298,TR);                               %Enthalpy departure of D [=] kJ/kmol
    ICP(5)= integral(CPE,298,TR);                               %Enthalpy departure of E [=] kJ/kmol
    ICP(6)= integral(CPI,298,TR);                               %Enthalpy departure of I [=] kJ/kmol
    

end

