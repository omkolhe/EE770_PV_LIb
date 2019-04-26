% Name  - Om Kolhe
% Assignment 3 -  EE770

function F = SDM_3eqn(x)

%% Values from the datasheet - Trina Solar (The TALLMAX Module) ? Module with 325 W.

%Constants
N_Cell = 72;
Vt = 25.7e-3; % Thermal voltage at temperature 25 degree C.

% Values at STC 
V_OC_STC = 45.9; % in Volts
I_SC_STC = 9.25; % in Amperes
I_MP_STC = 8.76; % in Amperes
V_MP_STC = 37.2; % in Volts

%Input x vector Describtion - 
% x(1) = Rs
% x(2) = Rsh
% x(3) = n

expo = exp((-V_OC_STC + V_MP_STC + (I_MP_STC*x(1)))/(x(3)*N_Cell*Vt));

frac = (-V_OC_STC + (I_SC_STC*x(2))+ (I_SC_STC*x(1)))/(x(3)*N_Cell*Vt*x(2));

expo2 = exp((-V_OC_STC + (I_SC_STC*x(1)))/(x(3)*N_Cell*Vt));

F(1,1) = I_MP_STC - I_SC_STC + (V_MP_STC + (I_MP_STC*x(1)) - (I_SC_STC*x(1)))/x(2) + (I_SC_STC - (V_OC_STC-I_SC_STC*x(1))/x(2))*expo;

F(2,1) = I_MP_STC + V_MP_STC*(((-frac*expo)- 1/x(2))/(1 + (frac*x(1)*expo) + x(1)/x(2)));

F(3,1) = 1/x(2) + ((-frac*expo2) - 1/x(2))/(1 + (frac*x(1)*expo2) + x(1)/x(2));

end

