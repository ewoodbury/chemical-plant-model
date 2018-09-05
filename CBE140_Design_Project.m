function [Condenser_temperature, Minimum_cost] = CBE140_Design_Project()
%Calculates the optimal operating temperature (in C) of condenser to minimize
%annual cost of producing ammonia using a process that follows the design specs.

X_effective_vector = linspace(.6,.954,500);
Cost_vector = zeros(1,500);
Temperature_vector = zeros(1,500);

%Setting up for loop to loop through entire range of X_effective
for i = 1:500
    Volume_cost = Volume_cost_from_X_effective(X_effective_vector(i));
    Cooling_cost = Cooling_cost_from_X_effective(X_effective_vector(i));
    Cost_vector(i) = Volume_cost + Cooling_cost;
    Temperature_vector(i) = T_from_X_effective(X_effective_vector(i));
end

%Finding the minimum cost and the temperature associated with minimum cost:
[Minimum_cost, index] = min(Cost_vector);
Condenser_temperature = Temperature_vector(index);

%Plotting total cost versus condenser temperature:
plot(Temperature_vector, Cost_vector);
title('Total Cost versus Condenser Temperature')
xlabel('Condenser Temperature (°Celcius)')
ylabel('Total Cost ($/year)')
axis([-83 60 2.1e8 2.7e8])
end

function [Volume_cost] = Volume_cost_from_X_effective(X_effective)
%Calculates the annual cost of operating a reactor from X_effective

%Finding Volume of reactor
V = Volume_from_X_effective(X_effective);

%Finding annual cost from given price $200000/ft^3
Volume_cost = 200000.*V;
end

function [V] = Volume_from_X_effective(X_effective)
%Gives volume of the reactor from the overall conversion and gas mole
%fractions in stream 5.

%Finding gas mole fractions in stream 7:
[yN2,~,~,yNH3] = y_from_X_effective(X_effective);

%Finding F1 and F2:
F1 = 2.69221e3./X_effective;
F2 = 5.*F1;
F7 = 4.*F1;

%Define moles of each N2, H2, NH3 going into rxr (stream 2):
nN2 = .249275.*F1 + yN2.*F7;
nH2 = 3.*nN2;
nNH3 = yNH3.*F7;


%Define the mol fractions of gases as function of single-pass conversion
%varies (X_s = single-pass conversion):
xN2 = @(X_s) (nN2 - X_s.*nN2)./F2;
xH2 = @(X_s) (nH2 - X_s.*nH2)./F2;
xNH3 = @(X_s) (nNH3 + 2.*X_s.*nN2)./F2;

%Define constants necessary for r_N2:
k1 = 1.2;
K = 4.67e-3;
P = 300; %bar
beta = .00140;

%Define r_N2:
r_N2 = @(X_s) -k1.*(xN2(X_s).*xH2(X_s).^3 - xNH3(X_s).^3./(K^2*P^2))./(xNH3(X_s).*xH2(X_s).^(1/2) + beta*xH2(X_s).^2).^(3/2);

%Define differential equation using PFR design equation dX/dV = -r/n_in
%(where n_in = nN2):
dVdX_s = @(X_s, V) nN2./-r_N2(X_s);

%Specifying initial condition and X_s_boundaries:
initial_V = 0;
X_s_boundaries = [0,.2];

%Using ode45 to solve diff eq:
[~, Volume_range] = ode45(dVdX_s, X_s_boundaries, initial_V);

%Volume of reactor is the V at X_2 = .2, which is the last value in the
%table:
V = Volume_range(end);
end

function [yN2, yH2, yA, yNH3] = y_from_X_effective(X_effective)
%Calculates vapor mole fractions in the condenser from X_effective
%final temp should be 279.6K (+/- 20K is fine)
F4 = X_effective.*2*.249275;
f = @(y) 1.6.*y.^2 + (-1.7003 + F4).*y + .19942;
yN2 = fsolve(f,.1);
F6 = (.9003 - F4) - 1.6.*yN2;
yH2 = 3.*yN2;
yA = .0029./F6;
yNH3 = ((.4*.249275) - F4 + 1.6.*yN2)./F6;
end

function [Cooling_cost] = Cooling_cost_from_X_effective(X_effective)
%Calculates the annual cooling cost from the X_effective and F3 (F3 is
%lbmol/hr)

%Finding F3:
F3 = F3_from_X_effective(X_effective);

%Finding total lbmol that flows through F3 per year:
F3_year = F3*24*365;

%Finding temperature that condenser is operating at:
T_condenser = T_from_X_effective(X_effective);

%Defining molar heat capacity of the gas stream out of reactor and reactor temperature:
Cp = 17.74; %(units = Btu/(lbmol*K)
T_reactor = (900 - 32)*(5/9); %have to convert from 900F to C

%Using given cooling cost equation to find BTU used annually:
Q = F3_year*Cp*(T_reactor - T_condenser);

%Calculating cost from given price $.0002/Btu
Cooling_cost = .0002*Q;
end

function [F3] = F3_from_X_effective(X_effective)
%Calculates the molar flow rate of stream 3 given X_effective

%Finding the vapor mole fractions in F5,F6,F7:
[yN2,~, yA, yNH3] = y_from_X_effective(X_effective);

%Finding F1,F7:
F1 = 2.69221e3./X_effective;
F7 = 4.*F1;

%Finding molar flow rates of each component in F2:
n2N2 = .249275.*F1 + yN2.*F7;
n2H2 = 3.*n2N2;
n2A = .0029*F1 + yA*F7;
n2NH3 = yNH3.*F7;

%Using 20% single-pass conversion to calculate component flow rates in F3:
n3N2 = .8*n2N2;
n3H2 = .8*n2H2;
n3A = n2A;
n3NH3 = n2NH3 + 2*.2*n2N2;

%Summing components to calculate total F3 flow rate (in lbmol/hr)
F3 = n3N2 + n3H2 + n3A + n3NH3;
end

function [T_condenser] = T_from_X_effective(X_effective)
%Calculates condenser temperature given X_effective

%Calculating vapor mole fraction of ammonia in condenser:
[~,~,~, yNH3] = y_from_X_effective(X_effective);

%Using vapor mole fraction to calculate condenser temperature:
T_condenser = T_from_yNH3(yNH3);
end

function [T] = T_from_yNH3(yNH3)
%This function takes the mole fraction of ammonia and gives the temperature
%of the condenser in Celsius.

% Establishing Antoine Constants:
A = 7.36050;
B = 926.132;
C = 240.17;
%Putting system pressure in mmHg:
P = 300 * 760 / 1.01325;
%Using rearranged Antoine equation to find T:
T = (B / (A - log10(yNH3.*P))) - C;
end

