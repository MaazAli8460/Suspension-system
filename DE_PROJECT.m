printf('Name of software house: FAST ISB \nTEAM MEMBERS: \n');
printf('Name          : Student Id  :\n');


choice=1;
while(choice!=0)
key=input('Enter any key to continue:');
m = 1800; % mass of the satellite in kg
k1 = 220; % spring constant in N/m for Air Coasters Shock
B1 = 100; % damping coefficient in N.s/m for Air Coasters Shock
x0 = 0.03; % displacement before loading in m
L1=B1./(2.*m);
% Analytical solution
W1 = sqrt(k1/m); % natural frequency in rad/s
c1 = x0; % constant of integration
c2 = 0.0068; % constant of integration
T1 = 0:0.1:100; % time interval for plotting
x = e.^(-L1.*T1).*(c1*cos((sqrt(W1.^2 -L1.^2)*T1)) + c2*sin((sqrt(W1.^2 -L1.^2)*T1))); % displacement of the shock
% Plotting the analytical solution
plot(T1, x);
xlabel('Time (s)');
ylabel('Displacement (m)');
title('Analytical Solution for Air Coasters Shock');
printf('\nThis is the Analytical solution for  Air Coasters shocks\n');
key=input('Enter any key to continue:');
% MATLAB solution
TRANGE =0:0.1:100; % time interval for solving ODE
v0=0;
cond = [x0 v0]; % initial conditions
[t, xa] = ode45(@(t,xa) [xa(2); (-B1/m)*xa(2) - (k1/m)*xa(1)], TRANGE, cond);

% Plotting the OCTAVE solution
figure;
xtempa=xa(:,1);
plot(t, xtempa);
xlabel('Time (s)');
ylabel('Displacement (m)');
title('OCTAVE Solution for Air Coasters Shock');
printf('\nThis is the OCTAVE solution for  Air Coasters shocks\n');
key=input('Enter any key to continue:');

k2 = 59; % spring constant of 2nd type of shocks "Cloud Raiders".
B2 = 1050; % damping coefficient of 2nd type of shocks "Cloud Raiders".
x0 = 0.03; % displacement before loading in m

% Analytical solution
L2 =B2./(2.*m);% lemda of 2nd type of shocks "Cloud Raiders".
W2 = sqrt(k2/m); % omega of 2nd type of shocks "Cloud Raiders".
c1 = 0.032;
c2 = -0.003;
T2 = 0:0.1:100; % time interval for plotting the graph
x = (e.^(-L2.*T2)).*((c1.*e.^(sqrt(-W2.^2 +L2.^2).*T2)) + (c2.*e.^(sqrt(-W2.^2 +L2.^2).*T2))); % displacement of the shock

% Plotting the analytical solution

plot(T2, x);
xlabel('Time in seconds (s)');
ylabel('Displacement in meters (m)');
title('Analytical Solution for Cloud Raiders Shock');
printf('\nThis is the analytical solution for Cloud Raiders shocks\n');

% OCTAVE solution
key=input('Enter any key to continue:');
TRANGE = 0:0.1:100; % time interval for solving ODE
v0=0;
cond = [x0 v0]; % initial conditions

[T, xb] = ode45(@(T,xb) [xb(2); (-B2/m)*xb(2) - (k2/m)*xb(1)], TRANGE, cond);

% Plotting the OCTAVE solution
figure;
xtempb=xb(:,1);
plot(T, xtempb);
xlabel('Time (s)');
ylabel('Displacement (m)');
title('OCTAVE Solution for Cloud Raiders Shock');
printf('\nThis is the Octave solution for Cloud Raiders shocks using ODE 45 function\n');

key=input('Enter any key to continue:');

k3 = 170; % spring constant of 2nd type of shocks "Cloud Raiders".
B3 = 1110; % damping coefficient of 2nd type of shocks "Cloud Raiders".
x0 = 0.03; % displacement before loading in m

% Analytical solution
L3 =B3./(2.*m);% lemda of 3nd type of shocks "EXTREME Shocks".
W3 = sqrt(k3/m); % omega of 3nd type of shocks "Extreme Shocks".
c1 = 0.2;
c2 = -0.17;
T3 = 0:0.1:100; % time interval for plotting the graph
x = (e.^(-L3.*T3)).*((c1.*e.^(sqrt(-W3.^2 +L3.^2).*T3)) + (c2.*e.^(sqrt(-W3.^2 +L3.^2).*T3))); % displacement of the shock

% Plotting the analytical solution

plot(T3, x);
xlabel('Time in seconds (s)');
ylabel('Displacement in meters (m)');
title('Analytical Solution for EXTREME Shock');
printf('\nThis is the analytical solution for EXTREME shocks\n');

% OCTAVE solution
key=input('Enter any key to continue:');
TRANGE = 0:0.1:100; % time interval for solving ODE
v0=0;
cond = [x0 v0]; % initial conditions

[T1, xc] = ode45(@(T1,xc) [xc(2); (-B3/m)*xc(2) - (k3/m)*xc(1)], TRANGE, cond);

% Plotting the OCTAVE solution
figure;
xtempc=xc(:,1);
plot(T1, xtempc);
xlabel('Time (s)');
ylabel('Displacement (m)');
title('GNU OCTAVE Solution for EXTREME Shock');
printf('\nThis is the Octave solution for EXTREME shocks using ODE 45 function\n');
key=input('Enter any key to continue:');
%this is the combined plot overlay of 3 types of shocks using ODE 45 function of GNU OCTAVE
printf('\nCombined Plot of OCTAVE GNU ODE-45 FUNCTION\n');
plot(t,xtempa, 'b');
xlabel('Time in seconds (s)');
ylabel('Displacement in meters (m)');

hold on;
plot(T,xtempb, 'r');

hold on;
plot(T1,xtempc, 'g');
legend('Air Coasters','Cloud Raiders','Extreme Shocks');
choice=0;
choice=input('If you want to QUIT press 0 else to continue press 1:');
endwhile;
