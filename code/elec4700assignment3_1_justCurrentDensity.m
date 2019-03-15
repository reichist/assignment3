clear
clc

% Variables
global C

C.m0 = 9.11e-31;         % rest mass of electron
C.mn = 0.26*C.m0;        % effective mass of electron
C.k = 1.381e-23;         % Boltzmann constant
C.q = 1.6023e-19;         % charge of electron

width = 200e-9;          % width of region
height = 100e-9;         % height of region
T = 300;                 % temperature
v_th = sqrt((C.k * T )/ C.mn); % thermal velocity
tau_mn = 0.2e-12;       % mean time between collisions
lambda = v_th * tau_mn; % mean free path
voltage_x = 0.1;        % applied voltage in x
voltage_y = 0;          % applied voltage in y
time_interval = 5e-15;  % time between steps in seconds

n = 1000;             % number of electrons simulated
steps = 1000;         % number of steps simulated

% Electron Simulation
% Vector Setup
electrons_x = rand(1, n)*width;
electrons_y = rand(1, n)*height;
electrons_vx = (v_th/sqrt(2)).*randn(1, n);
electrons_vy = (v_th/sqrt(2)).*randn(1, n);

new_electrons_x = zeros(1, n);
new_electrons_y = zeros(1, n);
old_z = 0;
old_temperature = 300;
old_current_density = 0;

% Scattering Setup
p_scat = 1-exp(-time_interval/tau_mn);

% Temperature Distribution Setup
temperature_matrix = zeros(11, 11);

% Electric Field Setup/Calculation
electric_field_x = voltage_x/width;
electric_field_y = voltage_y/height;
fprintf('We can rearrange the equation V = E/d to calculate E\n');
fprintf('Vx = %3.3d, d = %3.3d\n', voltage_x, width);
fprintf('Vy = %3.3d, d = %3.3d\n', voltage_y, height);
fprintf('The electric field on the electrons in the x direction is %3.3dV/m\n', electric_field_x);
fprintf('The electric field on the electrons in the y direction is %3.3dV/m\n\n', electric_field_y);

% Electric Field Force Calculation
force_x = electric_field_x*C.q;
force_y = electric_field_y*C.q;
fprintf('We can calculate the force on each electron with the equation F = E*q\n');
fprintf('Ex = %3.3d, q = %3.3d\n', electric_field_x, C.q);
fprintf('Ey = %3.3d, q = %3.3d\n', electric_field_y, C.q);
fprintf('The force on each electron in the x direction is %3.3dN\n', force_x);
fprintf('The force on each electron in the y direction is %3.3dN\n\n', force_y);

% Acceleration of each electron
acceleration_x = force_x/C.mn;
acceleration_y = force_y/C.mn;
fprintf('We can rearrange the equation F = m*a to calculate F\n');
fprintf('Fx = %3.3d, m = %3.3d\n', force_x, C.mn);
fprintf('Fy = %3.3d, m = %3.3d\n', force_y, C.mn);
fprintf('The acceleration on the electrons in the x direction is %3.3dm/s2\n', acceleration_x);
fprintf('The acceleration on the electrons in the y direction is %3.3dm/s2\n\n', acceleration_y);

% Electron Drift Current Density/Average Carrier Velocity
concentration = 1e15;
concentration_m2 = concentration/1e-4;

figure(2)
clf

for z = 1:steps
    % Check for random scattering
    a=rand(1, n);
    electrons_vx(a<p_scat) = (v_th/sqrt(2))*randn(1, length(electrons_vx(a<p_scat))); 
    electrons_vy(a<p_scat) = (v_th/sqrt(2))*randn(1, length(electrons_vx(a<p_scat))); 
    
    % Update velocity of each electron
    electrons_vx = electrons_vx + time_interval*acceleration_x;
    electrons_vy = electrons_vy + time_interval*acceleration_y;
    
    % New X&Y position calculations
    new_electrons_x = electrons_x + time_interval*electrons_vx;
    new_electrons_y = electrons_y + time_interval*electrons_vy;

    % Check for BCs
    index = new_electrons_x>width;
    new_electrons_x(index) = new_electrons_x(index) - width;
    electrons_x(index) = electrons_x(index) - width;
    
    index = new_electrons_x<0;
    new_electrons_x(index) = new_electrons_x(index) + width;
    electrons_x(index) = electrons_x(index) + width;
    
    index = new_electrons_y>height;
    electrons_vy(index) = -electrons_vy(index);
    
    index = new_electrons_y<0;
    electrons_vy(index) = -electrons_vy(index);
        
    % Electron Drift Current Density
    drift_velocity = mean(electrons_vx);
    current_density = C.q*concentration*drift_velocity;
    
    % Plotting Current Density
    plot([old_z z], [old_current_density current_density], 'r');
    title('Current Density');
    xlabel('Time Step'); ylabel('Current Density (A/cm2)');
    xlim([1 steps]);
    hold on;

    pause(0.01);
    
    % Update electron coordinates
    electrons_x = new_electrons_x;
    electrons_y = new_electrons_y;
    old_z = z;
    old_current_density = current_density;
    fprintf('Step #%3d (%4.1f Complete)\n', z, (z/steps)*100);
end