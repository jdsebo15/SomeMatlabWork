function CPW = heatCapW( THX )

%M file which calculates the water heat capacity in the heat exchanger

    CPW(1)= (103244 - 173.15*(THX) + (.2681)*(THX)^2)/1000/18;
    %Water Heat Capacity [=] kJ/(kg*K)

end

