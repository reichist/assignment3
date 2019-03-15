% Clear everything
clear all
clc

global C 
% Variables declaration
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
voltage = 0.1;          % applied voltage
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
total_temperature = 0;
old_total_temperature = 0;
old_current_density = 0;

% Scattering Setup
p_scat = 1-exp(-time_interval/tau_mn);

% Temperature Distribution Setup
temperature_matrix = zeros(11, 11);

% Box Setup
 % cell 1 is left side of box
 % cell 2 is right side of box
 % cell 3 is bottom side of box
 % cell 4 is top side of box
box1coords = [80e-9 120e-9 1e-9 40e-9];
box2coords = [80e-9 120e-9 60e-9 100e-9];

% Check for Box spawns
 % box 1
for z = 1:n
    index1 = electrons_x > box1coords(1);
    index2 = electrons_x < box1coords(2);
    index3 = bitand(index1, index2);
    index4 = electrons_y < box1coords(4);
    index5 = bitand(index3, index4);
    if any(index5) == true
        replaced_x = rand(1, n)*width;
        electrons_x(index5) = replaced_x(index5);
    else
        break
    end
end

 % box 2
for z = 1:n
    index1 = electrons_x > box2coords(1);
    index2 = electrons_x < box2coords(2);
    index3 = bitand(index1, index2);
    index4 = electrons_y > box1coords(3);
    index5 = bitand(index3, index4);
    if any(index5) == true
        replaced_x = rand(1, n)*width;
        electrons_x(index5) = replaced_x(index5);
    else
        break
    end
end

% FD Setup
nx = width*1e9;
ny = height*1e9;

G = sparse(nx*ny, nx*ny);
B = zeros(nx*ny, 1);
sigma = ones(nx, ny);
sigma_in = 0.01;
sigma_out = 1;

lBC = 0.1;
rBC = 0;
uBC = 0;
dBC = 0;

fn = @(i, j) j + (i-1)*ny;

sigma((box1coords(1)*1e9):(box1coords(2)*1e9), (box1coords(3)*1e9):(box1coords(4)*1e9)) = sigma_out;
sigma((box2coords(1)*1e9):(box2coords(2)*1e9), (box2coords(3)*1e9):(box2coords(4)*1e9)) = sigma_out;

% Map G Matrix
for z = 1:nx
    for p = 1:ny
        location = fn(z, p);
        location_xm = fn(z-1, p);
        location_xp = fn(z+1, p);
        location_ym = fn(z, p-1);
        location_yp = fn(z, p+1);
        if z == 1
            G(location, location) = 1;
            B(location) = lBC;
        elseif z == nx
            G(location, location) = 1;
            B(location) = rBC;
        elseif p == 1
            sxp = (sigma(z,p)+sigma(z+1,p))/1.0;
            sxm = (sigma(z,p)+sigma(z-1,p))/1.0;
            syp = (sigma(z,p)+sigma(z,p+1))/2.0;

            G(location, location) = -(sxp+sxm+syp);
            G(location, location_xp) = sxp;
            G(location, location_xm) = sxm;
            G(location, location_yp) = syp;
        elseif p == ny
            sxp = (sigma(z, p) + sigma(z+1, p))/1.0;
            sxm = (sigma(z, p) + sigma(z-1, p))/1.0;
            sym = (sigma(z, p) + sigma(z, p-1))/2.0;

            G(location, location) = -(sxp + sxm + sym);
            G(location, location_xp) = sxp;
            G(location, location_xm) = sxm;
            G(location, location_ym) = sym;
        else
            sxm = (sigma(z,p)+sigma(z-1,p))/2.0;
            sxp = (sigma(z,p)+sigma(z+1,p))/2.0;
            sym = (sigma(z,p)+sigma(z,p-1))/2.0;
            syp = (sigma(z,p)+sigma(z,p+1))/2.0;

            G(location, location) = -(sxm+sxp+sym+syp);
            G(location, location_xp) = sxp;
            G(location, location_xm) = sxm;
            G(location, location_yp) = syp;
            G(location, location_ym) = sym;
        end
    end
end

X = G\B;
voltage_distribution = zeros(nx, ny);

% Remap Voltage Matrix
for z = 1:nx
    for p = 1:ny
        location = fn(z, p);
        voltage_distribution(z, p) = X(location);
    end
end

voltage_distribution = voltage_distribution';

% Electric Field Setup/Calculation
dx = width/nx;
dy = height/ny;
[electric_field_x, electric_field_y] = gradient(voltage_distribution);
electric_field_x = -electric_field_x;
electric_field_y = -electric_field_y;
electric_field_x = electric_field_x'/dx;
electric_field_y = electric_field_y'/dy;

% Electric Field Force Calculation
force_x = electric_field_x*C.q;
force_y = electric_field_y*C.q;

% Acceleration of each electron
acceleration_x = force_x/C.mn;
acceleration_y = force_y/C.mn;

% Electron Drift Current Density/Average Carrier Velocity
concentration = 1e15; % 10^15 cm^-2
concentration_m2 = concentration/1e-4;

figure(1)
clf

% Simulation
for z = 1:steps
    % Check for random scattering
    a = rand(1, n);
    electrons_vx(a<p_scat) = (v_th/sqrt(2))*randn(1, length(electrons_vx(a<p_scat))); 
    electrons_vy(a<p_scat) = (v_th/sqrt(2))*randn(1, length(electrons_vx(a<p_scat))); 
    
    % Update X and Y velocities
    for x = 1:nx
        index1 = electrons_x > (x-1)*1e-9;
        index2 = electrons_x <= (x)*1e-9;
        index3 = bitand(index1, index2);          
        for y = 1:ny
            index4 = electrons_y > (y-1)*1e-9;
            index5 = electrons_y <= (y)*1e-9;
            index6 = bitand(index4, index5);       
            index7 = bitand(index3, index6);     
            electrons_vx(index7) = electrons_vx(index7) + time_interval*acceleration_x(x, y); 
            electrons_vy(index7) = electrons_vy(index7) + time_interval*acceleration_y(x, y);  
        end
    end
    
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
    
    % Check for box BCs
     % Box 1
      % left side
    index1 = electrons_y < box1coords(4);         % check for electron below the top of box 1
    index2 = electrons_y > box1coords(3);         % check for electron above the bottom of box 1
    index3 = electrons_x < box1coords(1);         % check for old electron location on the left side of the left side of box 1
    index4 = new_electrons_x > box1coords(1);     % check for new electron location on the right side of the left side of box 1
    index5 = bitand(index1, index2);              % compare for y coordinate -> box1coords(3) < electrons_y < box1coords(4)
    index6 = bitand(index3, index4);              % compare for x coordinate -> electrons_x < box1coords(1) < new_electrons_x
    index7 = bitand(index5, index6);              % compare for qualifying electrons 
    electrons_vx(index7) = -electrons_vx(index7);
    
      % right side
    index1 = electrons_y < box1coords(4);         % check for electron below the top of box 1
    index2 = electrons_y > box1coords(3);         % check for electron above the bottom of box 1
    index3 = electrons_x > box1coords(2);         % check for old electron location on the right side of the right side of box 1
    index4 = new_electrons_x < box1coords(2);     % check for new electron location on the left side of the right side of box 1
    index5 = bitand(index1, index2);              % compare for y coordinate -> box1coords(3) < electrons_y < box1coords(4)
    index6 = bitand(index3, index4);              % compare for x coordinate -> electrons_x > box1coords(2) > new_electrons_x
    index7 = bitand(index5, index6);              % compare for qualifying electrons 
    electrons_vx(index7) = -electrons_vx(index7);
    
      % top side
    index1 = electrons_x > box1coords(1);         % check for electron to the right of the left side of box 1
    index2 = electrons_x < box1coords(2);         % check for electron to the left of the right side of box 1
    index3 = electrons_y > box1coords(4);         % check for old electron location above the top side of box 1
    index4 = new_electrons_y < box1coords(4);     % check for new electron location below the top side of box 1
    index5 = bitand(index1, index2);              % compare for y coordinate -> box1coords(1) < electrons_x < box1coords(2)
    index6 = bitand(index3, index4);              % compare for x coordinate -> electrons_y > box1coords(4) > new_electrons_y
    index7 = bitand(index5, index6);              % compare for qualifying electrons 
    electrons_vy(index7) = -electrons_vy(index7);

     % Box 2
      % left side
    index1 = electrons_y < box2coords(4);         % check for electron below the top of box 2
    index2 = electrons_y > box2coords(3);         % check for electron above the bottom of box 2
    index3 = electrons_x < box2coords(1);         % check for old electron location on the left side of the left side of box 2
    index4 = new_electrons_x > box2coords(1);     % check for new electron location on the right side of the left side of box 2
    index5 = bitand(index1, index2);              % compare for y coordinate -> box2coords(3) < electrons_y < box2coords(4)
    index6 = bitand(index3, index4);              % compare for x coordinate -> electrons_x < box2coords(1) < new_electrons_x
    index7 = bitand(index5, index6);              % compare for qualifying electrons 
    electrons_vx(index7) = -electrons_vx(index7);
    
      % right side
    index1 = electrons_y < box2coords(4);         % check for electron below the top of box 2
    index2 = electrons_y > box2coords(3);         % check for electron above the bottom of box 2
    index3 = electrons_x > box2coords(2);         % check for old electron location on the right side of the right side of box 2
    index4 = new_electrons_x < box2coords(2);     % check for new electron location on the left side of the right side of box 2
    index5 = bitand(index1, index2);              % compare for y coordinate -> box2coords(3) < electrons_y < box2coords(4)
    index6 = bitand(index3, index4);              % compare for x coordinate -> electrons_x > box2coords(2) > new_electrons_x
    index7 = bitand(index5, index6);              % compare for qualifying electrons 
    electrons_vx(index7) = -electrons_vx(index7);
    
      % bottom side
    index1 = electrons_x > box2coords(1);         % check for electron to the right of the left side of box 2
    index2 = electrons_x < box2coords(2);         % check for electron to the left of the right side of box 2
    index3 = electrons_y < box2coords(3);         % check for old electron location below the bottom side of box 2
    index4 = new_electrons_y > box2coords(3);     % check for new electron location above the bottom side of box 2
    index5 = bitand(index1, index2);              % compare for y coordinate -> box2coords(1) < electrons_x < box2coords(2)
    index6 = bitand(index3, index4);              % compare for x coordinate -> electrons_y > box2coords(4) > new_electrons_y
    index7 = bitand(index5, index6);              % compare for qualifying electrons 
    electrons_vy(index7) = -electrons_vy(index7);   
        
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
    
    % Boxes Setup
    box1 = rectangle('Position', [80e-9 0     40e-9 40e-9]); % Creates Box 1 from (80nm, 0nm) to (120nm, 40nm)
    box2 = rectangle('Position', [80e-9 60e-9 40e-9 40e-9]); % Creates Box 2 from (80nm, 60nm) to (120nm, 100nm)
    
    title('Electron Modelling');
    xlabel('Width'); ylabel('Height');
    grid on;
    xlim([0 200e-9]);
    ylim([0 100e-9]);
    
    % Current Density Calculations
    drift_velocity = mean(electrons_vx);
    current_density = C.q*concentration*drift_velocity;
    
    % Plotting Current Density
    subplot(235)
    plot([old_z z], [old_current_density current_density], 'r');
    title('Current Density');
    hold on;
    xlabel('Time Step'); ylabel('Current Density');
    xlim([1 steps]); %ylim([0 100]);
   
    % Drift Velocity to Temperature
    V(1, :) = sqrt(electrons_vx(1, :).^2 + electrons_vy(1, :).^2);
    V_mean = mean(V.^2);
    temperature = V_mean*C.mn/C.k;
    total_temperature = total_temperature + temperature;

    % Plot Temperature Lines
    subplot(232)
    plot([old_z z], [old_temperature temperature], 'r');
    title('Temperature');
    hold on;
    xlabel('Time Step'); ylabel('Temperature (K)');
    xlim([1 steps]);

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
            temperature_matrix(5:6, 1:4) = 0;
            temperature_matrix(5:6, 7:10) = 0;
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
    
    pause(0.01);
    
    % Update electron coordinates
    electrons_x = new_electrons_x;
    electrons_y = new_electrons_y;
    old_z = z;
    old_temperature = temperature;
    old_current_density = current_density;
    fprintf('Step %4d/%4d\n', z, steps);
end

figure(2)
clf
surf(voltage_distributionT)
title('Voltage Distribution')
xlabel('nx'); ylabel('ny');
colorbar;

figure(3)
clf
subplot(211)
surf(electric_field_x)
title('Electric field - x component');
xlabel('nx'); ylabel('ny');
colorbar;

subplot(212)
surf(electric_field_y)
title('Electric field - y component');
xlabel('nx'); ylabel('ny');
colorbar;

figure(4)
clf
subplot(211)
surf(force_x)
title('Force in x direction');
xlabel('nx'); ylabel('ny');
colorbar; 

subplot(212)
surf(force_y)
title('Force in y direction');
xlabel('nx'); ylabel('ny');
colorbar;

figure(5)
clf
subplot(211)
surf(acceleration_x)
title('Acceleration in x direction');
xlabel('nx'); ylabel('ny');
colorbar;

subplot(212)
surf(acceleration_y)
title('Acceleration in y direction');
xlabel('nx'); ylabel('ny');
colorbar;

