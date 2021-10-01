function k = rateConstant(TR)

%M file which calculates the rate constants using the Arrhenius equation
%and constants provided in the problem statement

global k0 E R;

k0 = [ 70 , 9 , 1.5 , 5]*(10^3.35);          %Frequency factor, units to make rate [=] kmol/(kg*min)
E = [2.07 , 1.33 , 1.50 , 2.05]*(10^4.25);  %Units [=] J/mol
R = 8.314;                                   %Gas Constant [=] (J/(mol*K))
i= 1;

while i<=size(E,2)
    k(i) = k0(i)*exp(-E(i)/(R*TR));
    i=i+1;
end

end

