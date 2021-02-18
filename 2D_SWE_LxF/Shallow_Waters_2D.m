% Shallow Waters - TG2-2S
% Copyright (C) 2020 - Carlo De Falco - Marco Fois
%{
This script solves the shallow waters equations in Omega c R^2, 

dU     dF1     dF2
--  +  ---  +  ---  = R,
dt      dx      dy

where:

    [ h   ]       [     hu1            ]         [     hu2            ]
U = [ hu1 ]  F1 = [  hu1^2 + g/2(h^2)  ]   F2 =  [     hu1u2          ]
    [ hu2 ]       [     hu1u2          ]         [  hu2^2 + g/2(h^2)  ]

and the source R is given by

    [   0      ]
R = [ gh*dZ/dx ].
    [ gh*dZ/dy ]

We solve this problem with the 2-steps Taylor-Galerkin Method, using finite
elements of degree 0 (Q0) for the first step and degree 1 (Q1) for the
second one.

***************************************************************************

%}


% Cleaning
clear all
close all
clc




% *************************************************************************
% SPATIAL INFOS - for FEM (Q0,Q1)-elements
% *************************************************************************

%{ 
Domain extrema (lower-left corner and top-right corner):

                        o---------o (x1,y1)
                        |         |
                        |         |
                        |         |
                (x0,y0) o---------o

%}

% QUECEDO EXAMPLE - "A reappraisal of Taylor-Galerkin Algorithm for drying
% wetting areas in shallow water equation" - Quecedo, Pastor 
% Int. J. Numer. Meth. Fluids, 2002, 38:515-531.

x0 = -400;
y0 = 0;
x1 = 400;
y1 = 80;

Px = 60;
Py = 6;



mesh=square2mesh('s',x0,y0,x1,y1,Px,Py);
mesh = meshinfo(mesh);


% Lumped Mass matrix
L = Mass(mesh); 


%**************************************************************************
% TEMPORAL INFOS - II order Taylor exp. in time
%**************************************************************************

T = 1000;           % Final time
t = 0;              % Initial time
it = 1;             % No of iterations
dt = 2e-2;          % Initial time-step


%**************************************************************************
% COEFFICIENTS AND/OR USEFUL CONSTANTS/PARAMETERS
%**************************************************************************

g = 9.81;           % Gravity

%**************************************************************************
% INITIAL VALUES FOR U = [h, hu_x, hu_y], h and Z
%**************************************************************************

% Domain dofs
x = mesh.x';
y = mesh.y';

% Bottom-shape function Z: Omega c R^2 ---> R 
Z(1:mesh.NV,1) = zeros;
Z_x(1:mesh.NV,1) = zeros;
Z_y(1:mesh.NV,1) = zeros;

% Heigth function h: [0,T] x Omega c R^2 ---> R
h0 = 8 - sin(pi.*x/2/400);
% h0 = ones(size(x)); h0(x>=-400 & x<=-300) = 1.5;

% Velocity Field: U = (hu), where: u: [0,T] x Omega c R^2 ---> R^2
U0x(1:mesh.NV,1) =  0;
U0y(1:mesh.NV,1) =  0;

% Storing initial values
h = h0;
Ux = U0x;
Uy = U0y;

hsave = h0;
Uxsave = U0x;
Uysave = U0y;

%**************************************************************************
% CODE STARTS HERE
%**************************************************************************

while t < T
    
    t = t + dt;
    
    % First call to flux
    [Fhx, FU1x,FU2x] = F_x(h,Ux,Uy,g);
    [Fhy, FU1y,FU2y] = F_y(h,Ux,Uy,g);
    
    % Main call to Source function
    [~,SUx,SUy] = source(h,Z_x,Z_y,g);
    
    % Mean of h,U and fluxes on each element 
    for ip = 1:mesh.NP
            
        for vtx = 1:4
        
            V(vtx) = mesh.polygon(ip).vertices(vtx);
            
            
            % Mean of h, Ux, Uy
            h_mean(ip,1) = sum(h(V))/4;
            
            Ux_mean(ip,1) = sum(Ux(V))/4;
            Uy_mean(ip,1) = sum(Uy(V))/4;
            
           
            % Mean of SUx and SUy
            SUx_mean(ip,1) = sum(SUx(V))/4;
            SUy_mean(ip,1) = sum(SUy(V))/4;
            
            
            % we'll use these quantities later
            Z_x_mean(ip,1) = sum(Z_x(V))/4;
            Z_y_mean(ip,1) = sum(Z_y(V))/4;
            
        end
    end

%**************************************************************************
% FIRST STEP (TG2-2S): Q0 Finite element methods
%**************************************************************************

% Q0-FEM for flux and source 
[dh,dUx,dUy] = divF_q0(mesh,Fhx,FU1x,FU2x,Fhy,FU1y,FU2y);

    
    %{
        First step equations.
        Recall that:
    
        (*) Fhx ==> h*u_x                 (**) Fhy ==> h*u_y
        (*) FU1x ==> h*(u_x)^2 + g*h^2/2  (**) FU1y ==> h*(u_x)*(u_y)
        (*) FU2x ==> h*(u_x)*(u_y)        (**) FU2y ==> h*(u_y)^2 + g*h^2/2
        
        
  In this example we have no source!

%}
% 
%     v1 = Ux./h;
%     v2 = Uy./h;
%     v1(h<=0)=0;
%     v2(h<=0)=0;
%     Vel = [v1,v2];
%     
%     
%    
%   
%   DT = min(min(mesh.area)./(norm(vecnorm(Vel, 2 , 2),inf) + norm(sqrt(g*h),inf)));
   
        
    h_half = h_mean + (dt/2)*(dh)./(mesh.area(:));
    Ux_half = Ux_mean + (dt/2)*SUx_mean + (dt/2)*(dUx)./(mesh.area(:));
    Uy_half = Uy_mean + (dt/2)*SUy_mean + (dt/2)*(dUy)./(mesh.area(:));
    
   
%    surf(reshape([mesh.polygon(:).xb],Px,Py),reshape([mesh.polygon(:).yb],Px,Py),reshape(h_half,Px,Py));
%    axis off
%     drawnow


%**************************************************************************
% SECOND STEP (TG2-2S): Q1 Finite element methods
%**************************************************************************
   
% Second call to flux and source (with intermediate values)
[Fhx_half, FU1x_half,FU2x_half] = F_x(h_half,Ux_half,Uy_half,g);
[Fhy_half, FU1y_half,FU2y_half] = F_y(h_half,Ux_half,Uy_half,g); 
[~, SU1_half,SU2_half] = source(h_half,Z_x_mean,Z_y_mean,g);

[dxh,dxUx,dxUy,dyh,dyUx,dyUy] = dU(mesh,h,Ux,Uy,dt);

    
Fh_x = Fhx_half - dxh; Fh_y = Fhy_half - dyh;
FU1_x = FU1x_half - dxUx; FU1_y = FU1y_half - dyUx;
FU2_x = FU2x_half - dxUy; FU2_y = FU2y_half - dyUy;
    
    
 % Q1-FEM for flux and source (recall that R = 0!) 
[dh_f,dUx_f,dUy_f] = divF_q1(mesh,Fh_x,FU1_x,FU2_x,Fh_y,...
                        FU1_y,FU2_y);
                
h = h + L\(dt*dh_f);
Ux = Ux + L\(dt*dUx_f);
Uy = Uy + L\(dt*dUy_f);
        
%{ 
***************************************************************************
        BOUNDARY CONDITIONS: for this example we have
***************************************************************************    
                        Uy = 0, dUx/dy = 0, dh/dy = 0
     o----------------------------------------------------------------o
     |                                                                |
     |<-----            Ux = 0, dUy/dx = 0, dUx/dy = 0          ----->|
     |                                                                |
     o----------------------------------------------------------------o
                        Uy = 0, dUx/dy = 0, dh/dy = 0
     
     %}
     
     % Vector of the boundary nodes
     BN = [];
     
     for vt = 1:mesh.NV
         if mesh.vertex(vt).marker == 1    % In this way we fix only the vtx on the boundary
             BN = [BN vt];
         end
     end
          
for kk = BN
          
    if (abs(mesh.vertex(kk).x + 400)) < 10*eps(400) || (abs(mesh.vertex(kk).x - 400)) < 10*eps(400)
             % Now we find the left and right boundary
             Ux(kk) = 0;
    end

    if (abs(mesh.vertex(kk).y - 0)) < 10*eps(400) || (abs(mesh.vertex(kk).y -80)) < 10*eps(400)
             % Now we find the lower and upper boundary
             Uy(kk) = 0;
             
    end
         
end
     
   
%**************************************************************************
% OPTIMAL TIME STEP CHOICE
%**************************************************************************

    v1 = Ux./h;
    v2 = Uy./h;
    v1(h<=0)=0;
    v2(h<=0)=0;
    Vel = [v1,v2];
    alpha = 0.025;  
    
   
   dtold = dt;
   dt = alpha*min(min(mesh.area)./(norm(vecnorm(Vel, 2 , 2),inf) + norm(sqrt(g*h),inf)));
   %dt = 0.05*dt;
    dt = min(dt,2*dtold);

    
    it = it + 1;
    
     % Update main variables
     hsave = h;  
     Uxsave = Ux;
     Uysave = Uy;
     

figure(1)
surf(reshape(mesh.x,Px+1,Py+1),reshape(mesh.y,Px+1,Py+1),reshape(h,Px+1,Py+1));
title(['time [s] = ' num2str(t+dt)])
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]')
axis off
view(0,25)
drawnow  
pause(.001)


figure(2)
quiver(reshape(mesh.x,Px+1,Py+1),reshape(mesh.y,Px+1,Py+1),reshape(Ux,Px+1,Py+1),reshape(Uy,Px+1,Py+1));

     end
     
     


      
     
     
      





























