function CP = heatCap(TR)

%M file which calculates the per-species heat capacities inside the reactor

    CP(1)= (5.5 + .228*TR - (1.57*10^(-4))*TR^2);         %Heat Capacity of A [=] kJ/(kmol*K)
    CP(2)= (23.03 + .18*TR - (1.04*10^(-4))*TR^2);        %Heat Capacity of B [=] kJ/(kmol*K)
    CP(3)= (4.04 + .21*TR - (1.21*10^(-4))*TR^2);         %Heat Capacity of C [=] kJ/(kmol*K)
    CP(4)= (4.40 + .192*TR - (1.35*10^(-4))*TR^2);        %Heat Capacity of D [=] kJ/(kmol*K)
    CP(5)= (7.32 + .20*TR - (1.25*10^(-4))*TR^2);         %Heat Capacity of E [=] kJ/(kmol*K)
    CP(6)= (30.25 - (8.79*10^(-3)*TR) + (1.87*10^(-5)*TR^2) - (7.6*10^(-9)*TR^3));        %Heat Capacity of I [=] kJ/(kmol*K)
        
end

