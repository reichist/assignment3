% Clear everything
clear all
clc

% Variables
width = 200e-9;          % width of region
height = 100e-9;         % height of region
voltage = 0.1;           % applied voltage

nx = width*1e9;
ny = height*1e9;

sigma = ones(nx, ny);
sigma_in = .01; 
sigma_out = 1;

% Boxes Setup
box1coords = [80 120  1  40];
box2coords = [80 120 60 100];

G = sparse(nx*ny, nx*ny);
B = zeros(nx*ny, 1);

fn = @(i, j) j + (i-1)*ny;

sigma(box1coords(1):box1coords(2), box1coords(3):box1coords(4)) = sigma_in;
sigma(box2coords(1):box2coords(2), box2coords(3):box2coords(4)) = sigma_in;

for z = 1:nx
    for p = 1:ny
        n = fn(z, p);
        nxm = fn(z-1, p);
        nxp = fn(z+1, p);
        nym = fn(z, p-1);
        nyp = fn(z, p+1);
        if z == 1
            G(n, n) = 1;
            B(n) = 0.1;
        elseif z == nx
            G(n, n) = 1;
            B(n) = 0;
        elseif p == 1
            sxp = (sigma(z,p)+sigma(z+1,p))/1.0;
            sxm = (sigma(z,p)+sigma(z-1,p))/1.0;
            syp = (sigma(z,p)+sigma(z,p+1))/2.0;
            
            G(n, n) = -(sxp+sxm+syp);
            G(n, nxp) = sxp;
            G(n, nxm) = sxm;
            G(n, nyp) = syp;
        elseif p == ny
            sxp = (sigma(z, p) + sigma(z+1, p))/1.0;
            sxm = (sigma(z, p) + sigma(z-1, p))/1.0;
            sym = (sigma(z, p) + sigma(z, p-1))/2.0;
            
            G(n, n) = -(sxp + sxm + sym);
            G(n, nxp) = sxp;
            G(n, nxm) = sxm;
            G(n, nym) = sym;
        else
            sxm = (sigma(z,p)+sigma(z-1,p))/2.0;
            sxp = (sigma(z,p)+sigma(z+1,p))/2.0;
            sym = (sigma(z,p)+sigma(z,p-1))/2.0;
            syp = (sigma(z,p)+sigma(z,p+1))/2.0;
            
            G(n, n) = -(sxm+sxp+sym+syp);
            G(n, nxp) = sxp;
            G(n, nxm) = sxm;
            G(n, nyp) = syp;
            G(n, nym) = sym;
        end
    end
end

X = G\B;

voltage_distribution = zeros(nx, ny);
for z = 1:nx
    for p = 1:ny
        n = fn(z, p);
        voltage_distribution(z, p) = X(n);
    end
end

% Plot V(x,y)
figure(1)
clf
surf(voltage_distribution);
title('Voltage Distribution');
xlabel('Height (x10nm)'); ylabel('Width (x10nm)'); zlabel('Voltage (V)');
view(135, 45)
colorbar;

% Electric Field Setup/Calculation
[electric_field_x, electric_field_y] = gradient(voltage_distribution);

% Plot electric field vector plot
figure(2)
clf
quiver(electric_field_x.', electric_field_y.');
title('Electric Field');
xlabel('Width (nm)'); ylabel('Height (nm)');
grid on; 