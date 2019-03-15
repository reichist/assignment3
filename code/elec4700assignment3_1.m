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
fprintf('We can calculate E using the equation E = V/d\n');
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
fprintf('We can calculate accelerationusing the equation a = F/m\n');
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
        
    % Plot all electrons
    subplot(231)
    plot([electrons_x(1) new_electrons_x(1)], [electrons_y(1) new_electrons_y(1)], 'b');
    hold on;
    plot([electrons_x(2) new_electrons_x(2)], [electrons_y(2) new_electrons_y(2)], 'g');
    hold on;
    plot([electrons_x(3) new_electrons_x(3)], [electrons_y(3) new_electrons_y(3)], 'r');
    hold on;
    plot([electrons_x(4) new_electrons_x(4)], [electrons_y(4) new_electrons_y(4)], 'c');
    hold on;
    plot([electrons_x(5) new_electrons_x(5)], [electrons_y(5) new_electrons_y(5)], 'm');
    hold on;
    plot([electrons_x(6) new_electrons_x(6)], [electrons_y(6) new_electrons_y(6)], 'k');
    hold on;
    
    title('Electron Modelling');
    xlabel('Width'); ylabel('Height');
    grid on;
    xlim([0 200e-9]);
    ylim([0 100e-9]);
   
    % Drift Velocity to Temperature
    V(1, :) = sqrt(electrons_vx(1, :).^2 + electrons_vy(1, :).^2);
    V_mean = mean(V.^2);
    temperature = V_mean*C.mn/C.k;

    % Plot Temperature Lines
    subplot(232)
    plot([old_z z], [old_temperature temperature], 'r');
    title('Temperature');
    xlabel('Time Step'); ylabel('Temperature (K)');
    xlim([1 steps]);
    hold on;
    
    % Electron Density Map
    if rem(z, 10) == 0
        subplot(233)
        hist3([electrons_x', electrons_y'], 'CDataMode', 'auto', 'FaceColor', 'interp')
        colormap('default');
        colorbar;
        xlabel('Width'); ylabel('Height');
        title('Electron Density');
        view(2);
    end
    
    % Temperature Distribution Map
    if rem(z, 10) == 0
        for y = 1:10
            ymax = y*10;
            ymin = ymax-10;
            for x = 1:10
                xmax = x*20;
                xmin = xmax-20;
                index1 = electrons_x > (xmin*1e-9);
                index2 = electrons_x < (xmax*1e-9);
                index3 = electrons_y > (ymin*1e-9);
                index4 = electrons_y < (ymax*1e-9);
                index5 = bitand(index1, index2);
                index6 = bitand(index3, index4);
                index7 = bitand(index5, index6);
                velocity = sqrt((electrons_vx(index7).^2) + (electrons_vy(index7).^2));
                v_mean = mean(velocity.*velocity);
                temperature_value = (((v_mean)*(C.mn))/(C.k));
                temperature_matrix(x, y) = temperature_value;
            end       
        end
    end

     % Plotting Temperature Distribution
     if rem(z, 10) == 0
        subplot(234)
        surf(transpose(temperature_matrix));
        xlabel('Width (a.u.)'); ylabel('Height(a.u.)');
        title('Temperature Distribution');
        colormap('default');
        colorbar;
        view(2);
     end
        
    % Electron Drift Current Density
    drift_velocity = mean(electrons_vx);
    current_density = C.q*concentration*drift_velocity;
    
    % Plotting Current Density
    subplot(235)
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
    old_temperature = temperature;
    old_current_density = current_density;
    fprintf('Step #%4d/%4d\n', z, steps);
end