clear all; clc;

%% Problem Parameters
nx = 10; % Number of grid points in x-direction
ny = 10; % Number of grid points in y-direction
L = 1; % Length of the plate in x-direction
H = 1; % Height of the plate in y-direction
dx = L/(nx-1); % Grid spacing in x-direction
dy = H/(ny-1); % Grid spacing in y-direction
x = linspace(0,L,nx); % Grid points in x-direction
y = linspace(0,H,ny); % Grid points in y-direction
D = 1; % Diffusivity coefficient
U = 1; % Velocity of the fluid
V = 0; % Velocity of the fluid
Pe = U*L/D; % Peclet number

%% Initialization
T = zeros(nx,ny); % Temperature matrix
T(:,1) = 1; % Boundary condition at y=0
T(1,:) = 0; % Boundary condition at x=0
T(nx,:) = 0; % Boundary condition at x=L

%% Solution
figure;
v = VideoWriter('2D_convection_diffusion.avi');
open(v);
for t = 1:50
    for i = 2:nx-1
        for j = 2:ny-1
            % Convection term
            c1 = U*(T(i,j)-T(i-1,j))/dx;
            c2 = V*(T(i,j)-T(i,j-1))/dy;
            convective_term = Pe*(c1 + c2);

            % Diffusion term
            d1 = (T(i+1,j)-2*T(i,j)+T(i-1,j))/dx^2;
            d2 = (T(i,j+1)-2*T(i,j)+T(i,j-1))/dy^2;
            diffusive_term = D*(d1 + d2);

            % Update temperature
            T(i,j) = T(i,j) + (convective_term - diffusive_term);
        end
    end

    %% Plotting
    [X,Y] = meshgrid(x,y);
    surf(X,Y,T);
    title('2D Convection Diffusion on a Flat Plate');
    xlabel('x'); ylabel('y'); zlabel('Temperature');
    
    % Add the frame to the video
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);
