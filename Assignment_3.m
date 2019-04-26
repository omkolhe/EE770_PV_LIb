% Name  - Om Kolhe
% Assignment 3 -  EE770
%% Values from the datasheet - Trina Solar (The TALLMAX Module)  Module with 325 W.

%Constants
N_Cell = 72;
Vt = 25.7e-3; % Thermal voltage at temperature 25 degree C.

% Values at STC 
V_OC_STC = 45.9; % in Volts
I_SC_STC = 9.25; % in Amperes
I_MP_STC = 8.76; % in Amperes
V_MP_STC = 37.2; % in Volts

%% Initial Guesses 

n_initial = (2*V_MP_STC - V_OC_STC)/(N_Cell*Vt*(log((I_SC_STC-I_MP_STC)/I_SC_STC) + (I_MP_STC/(I_SC_STC-I_MP_STC))));
Rs_initial = (V_MP_STC/I_MP_STC) - ((2*V_MP_STC - V_OC_STC)/(I_SC_STC-I_MP_STC))/(log((I_SC_STC-I_MP_STC)/I_SC_STC) + (I_MP_STC/(I_SC_STC-I_MP_STC))); 
Rsh_initial = sqrt(Rs_initial/((I_SC_STC/(n_initial*N_Cell*Vt))*(exp((Rs_initial*I_SC_STC - V_OC_STC)/(n_initial*N_Cell*Vt)))));


%% Solving 3 non linear equations 

initial_guess = zeros(1,3);
initial_guess(1) = Rs_initial;
initial_guess(2) = Rsh_initial;
initial_guess(3) = n_initial;

SDMeqn = @SDM_3eqn;

options=optimset('TolFun',1e-6,'MaxIter',100000,'MaxFunEvals',100000);

x = fsolve(SDMeqn,initial_guess,options);

%Solutions of the above equations
Rs_STC = x(1);
Rsh_STC = x(2);
n_STC = x(3);

%Calculating Isat and Iph from the above results
Isat_STC = (I_SC_STC - (V_OC_STC - I_SC_STC*Rs_STC)/Rsh_STC)*(exp(-V_OC_STC/(n_STC*N_Cell*Vt)));

Iph_STC = Isat_STC*exp(V_OC_STC/(n_STC*N_Cell*Vt)) + V_OC_STC/Rsh_STC;

%% Plotting the IV curve at STC 

figure();
f = @(V,I) I - Iph_STC + Isat_STC*(exp((V+I*Rs_STC)/(n_STC*N_Cell*Vt))-1) + ((V+I*Rs_STC)/Rsh_STC);
ezplot(f, [0 50 0 10]);
title('IV Curve of the Module at STC');
ylabel('Current (in A)');
xlabel('Voltage (in V)');

%% Calculating the five parameters at different values of irradiance and temperatures

%Constants
Temp_STC = 25+273;  % in Kelvin
Temp_NOTC = 44+273; % in Kelvin
G_STC = 1000; % W per m2
G_NOTC = 800; % W per m2
I_SC_NOTC = 7.47; % in Amperes
alpha = 5e-4;
q = 1.60217662e-19;
E_G = 1.1; % for silicon 
K_B = 1.38064852e-23;

K = log((I_SC_NOTC - alpha*(Temp_NOTC-Temp_STC))/I_SC_STC)/(log(G_NOTC/G_STC));

%Initializing variables
I_SC = zeros(10,13);
Rs = zeros(10,13);
Rsh = zeros(10,13);
Iph = zeros(10,13);
Isat = zeros(10,13);
Voc = zeros(10,13);
Vmp = zeros(10,13);
Imp = zeros(10,13);
Pmp = zeros(10,13);
FF = zeros(10,13);

figure();

%Calculating the parameters at different values of irradiance and temperatures
for i = 1:10 
    G = 400 + (i-1)*100;
    
    for j = 1:13
        temp = 273 + 20 + (j-1)*5;
        I_SC(i,j) = ((G/G_STC)^K)*(I_SC_STC + alpha*(temp-Temp_STC));
        Rs(i,j) = Rs_STC;
        Rsh(i,j) = (G_STC/G)*Rsh_STC;
        Iph(i,j) = I_SC(i,j)*(1 + (Rs(i,j)/Rsh(i,j)));
        Isat(i,j) = Isat_STC*(temp/Temp_STC)^3*exp(((q*E_G)/(n_STC*K_B))*(1/Temp_STC - 1/temp));
        
        %Calculating Voc
        fun1 = @(V) Iph(i,j) - Isat(i,j)*(exp(V/(n_STC*N_Cell*Vt))-1) - V/Rsh(i,j);
        Voc(i,j) = fzero(fun1,40);
        
        %Finding the maximum power point in the IV curve
        Vt = K_B*temp/q; % Considering the variation of thermal voltage with temperature
        f1 = @(V,I) I-Iph(i,j)+Isat(i,j)*(exp((V+I*Rs(i,j))/(n_STC*N_Cell*Vt))-1)+((V+I*Rs(i,j))/Rsh(i,j)) ;
        IV = ezplot(f1, [0 60 0 15]);
        tmp = get(IV,'contourMatrix');
        V = tmp(1,3:end);
        I = tmp(2,3:end);
        P = V.*I;
        maximum = max(P);
        max_index = find(P==maximum);
        Vmp(i,j) = V(max_index);
        Imp(i,j) = I(max_index);
        Pmp(i,j) = maximum;
        FF(i,j) = Pmp(i,j)/(Voc(i,j)*I_SC(i,j));
    end
end 

%% IV Curves for different irradiances and temperatures

for temp = [2 5 7 9 11 13]
    temperature = 273 + 20 + (temp-1)*5;
    for i = 1:10
        Iph_1 = Iph(i,temp);
        Isat_1 = Isat(i,temp);
        Rs_1 = Rs(i,temp);
        Rsh_1 = Rsh(i,temp);
        Vt = K_B*temperature/q; % Considering the variation of thermal voltage with temperature
        f = @(V,I) I - Iph_1 + Isat_1*(exp((V+I*Rs_1)/(n_STC*N_Cell*Vt))-1) + ((V+I*Rs_1)/Rsh_1);
        ezplot(f, [0 50 0 15]);
        hold on;
    end
    title(['IV Curve for different irradiance at ' num2str(temperature-273) ' Degree Celsius']);
    ylabel('Current (in A)');
    xlabel('Voltage (in V)');
    %legend({'400','500','600','700','800','900','1000','1100','1200','1300'});
    figure();
    
end 

%% Plotting contour plots for different irradiances and temperatures
temper = 20:5:80;
irr = 400:100:1300;

figure();
surf(temper,irr,I_SC);
view(2);
title('Contour Plot of Short Circuit Current');
ylabel('Irradiance in W per m2');
xlabel('Temperature in degree Celsius');

figure();
surf(temper,irr,Voc);
view(2);
title('Contour Plot of Open Circuit Voltage');
ylabel('Irradiance in W per m2' );
xlabel('Temperature in degree Celsius');

figure();
surf(temper,irr,Pmp);
view(2);
title('Contour Plot of Maximum Power');
ylabel('Irradiance in W per m2' );
xlabel('Temperature in degree Celsius');

figure();
surf(temper,irr,FF);
view(2);
title('Contour Plot of Fill Factor');
ylabel('Irradiance in W per m2' );
xlabel('Temperature in degree Celsius');

figure();
surf(temper,irr,Rsh);
view(2);
title('Contour Plot of Shunt Resistance');
ylabel('Irradiance in W per m2');
xlabel('Temperature in degree Celsius');

figure();
surf(temper,irr,Iph);
view(2);
title('Contour Plot of Iph');
ylabel('Irradiance in W per m2');
xlabel('Temperature in degree Celsius');

figure();
surf(temper,irr,Isat);
view(2);
title('Contour Plot of Isat');
ylabel('Irradiance in W per m2');
xlabel('Temperature in degree Celsius');
