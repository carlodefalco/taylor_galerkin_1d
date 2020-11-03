% Peraire 1986 - Quecedo 2002 - Shallow Waters - TG2-2S
clear all
clc
close all
%{
This script solves the shallow waters equations in Omega c R^2, 

dU     dF1     dF2
--  +  ---  +  ---  = R,
dt      dx      dy

where:

    [ h   ]       [     hu1            ]         [     hu2            ]
U = [ hu1 ]  F1 = [hu1^2 + g/2(h^2-H^2)]   F2 =  [     hu1u2          ]
    [ hu2 ]       [     hu1u2          ]         [hu2^2 + g/2(h^2-H^2)]

and the source R is given by

    [   0      ]
R = [ gh*dZ/dx ].
    [ gh*dZ/dy ]

We solve this problem with the 2-steps Taylor-Galerkin Method, using finite
elements of degree 0 (Q0) for the first step and degree 1 (Q1) for the
second one.

***************************************************************************

%}


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
% Int. J. Numer. Meth. Fluids, 2002, 38:515-531
x0 = -400;
y0 = 0;
x1 = 400;
y1 = 80;

Nx = 60;
Ny = 6;

mesh=square2mesh('s',x0,y0,x1,y1,Nx,Ny);
mesh = meshinfo(mesh);

% Global Mass matrix 
M = Mass(mesh);
% Lumped Mass matrix
L = diag(M*ones(mesh.NV,1));



%**************************************************************************
% TEMPORAL INFOS - II order Taylor exp. in time
%**************************************************************************

T = 1500;           % Final time
t = 0;              % Initial time
it = 1;             % No of iterations
dt = 1e-2;          % Initial time-step


nsave = round(T/5);
tsave = linspace(0,T,nsave);

%**************************************************************************
% COEFFICIENTS AND/OR USEFUL CONSTANTS/PARAMETERS
%**************************************************************************

g = 9.81;           % Gravity

%**************************************************************************
% INITIAL VALUES FOR U = [h, hu_x, hu_y], h and Z
%**************************************************************************

% Domain dofs
x = mesh.xv';
y = mesh.yv';

% Bottom-shape function Z: Omega c R^2 ---> R 
Z(1:mesh.NV,1) = zeros;
Z_x(1:mesh.NV,1) = zeros;
Z_y(1:mesh.NV,1) = zeros;



% Heigth function h: [0,T] x Omega c R^2 ---> R
h0 = 8 - sin(pi.*x/2/400); 



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

for is = 2:nsave
     while (t(it) < tsave (is))

    t(it+1) = t(it) + dt;
    
    % This check that T is not exceeded
    if (t(it+1) > tsave (is))
        t(it+1) = tsave (is);
        dt = t(it+1) - t(it);
    end
    
    fprintf ("iteration No = %d; time progress = %g, dt = %g\n", it, t(it+1), dt)
    
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
            h_mean(ip,1) = sum(h(V))/mesh.area(ip);
            Ux_mean(ip,1) = sum(Ux(V))/mesh.area(ip);
            Uy_mean(ip,1) = sum(Uy(V))/mesh.area(ip);
            
           
            % Mean of SUx and SUy
            SUx_mean(ip,1) = sum(SUx(V))/mesh.area(ip);
            SUy_mean(ip,1) = sum(SUy(V))/mesh.area(ip);
            
            
            % we'll use these quantities later
            Z_x_mean(ip,1) = sum(Z_x(V))/mesh.area(ip);
            Z_y_mean(ip,1) = sum(Z_y(V))/mesh.area(ip);
            
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
        
        
     %}
    
        
    h_half = h_mean + (dt/2)*(dh)./(mesh.area(:));
    Ux_half = Ux_mean + (dt/2)*SUx_mean + (dt/2)*(dUx)./(mesh.area(:));
    Uy_half = Uy_mean + (dt/2)*SUy_mean + (dt/2)*(dUy)./(mesh.area(:));
    
    
    
    
    
%     
%    surf(reshape([mesh.polygon(:).xb],Nx,Ny),reshape([mesh.polygon(:).yb],Nx,Ny),reshape(dUx,Nx,Ny));
%     drawnow


%**************************************************************************
% SECOND STEP (TG2-2S): Q1 Finite element methods
%**************************************************************************
   
   % Second call to flux and source (with intermediate values)
    [Fhx_half, FU1x_half,FU2x_half] = F_x(h_half,Ux_half,Uy_half,g);
    [Fhy_half, FU1y_half,FU2y_half] = F_y(h_half,Ux_half,Uy_half,g); 
    [~, SU1_half,SU2_half] = source(h_half,Z_x_mean,Z_y_mean,g);
    
    
    % Q1-FEM for flux and source 
    [dh_f,dUx_f,dUy_f] = divF_q1(mesh,Fhx_half,FU1x_half,FU2x_half,Fhy_half,...
                        FU1y_half,FU2y_half);
                    
                    
         
    %{
        BOUNDARY CONDITIONS: for this example we have
     
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
     % Vector of the internal nodes
     IN = [];
     
     for vt = 1:mesh.NV
         if mesh.vertex(vt).marker == 1    % In this way we fix only the vtx on the boundary
             BN = [BN vt];
         else
             IN = [IN vt];
         end
     end
     

     
     for kk = BN
           
         if mesh.vertex(kk).x == -400
             % Now we find the left boundary
             Ux(kk) = 0;
         elseif mesh.vertex(kk).x == 400
             % Now we find the right boundary
             Ux(kk) = 0;
         elseif mesh.vertex(kk).y == 0
             % Now we find the lower boundary
             Uy(kk) = 0;
         elseif mesh.vertex(kk).y == 80
             % Now we find the upper boundary
             Uy(kk) = 0;
             
         end
         
     end
     
     %%%%%%%%%%%%% non ho riportato le sorgenti perché le ho su un altro m-file, 
     %%%%%%%%%%%%% ma nel nostro esempio sono zero (Z = 0) e
     %%%%%%%%%%%%% non portano contributo
      h = h + L\(dt*dh_f);
     Ux(IN) = Ux(IN) + L(IN,IN)\(dt*dUx_f(IN)); 
     Uy(IN) = Uy(IN) + L(IN,IN)\(dt*dUy_f(IN));
     
     
    
%**************************************************************************
% OPTIMAL TIME STEP CHOICE
%**************************************************************************

    v1 = Ux./h;
    v2 = Uy./h;
    v1(h<=0)=0;
    v2(h<=0)=0;
    Vel = [v1;v2];
    alpha = .9/sqrt(3);
    
    dtold = dt;
    dt = alpha*min(min(mesh.area)./(norm(Vel, inf) + norm(sqrt(g*h),inf)));
    dt = min(dt,2*dtold);
    
    it = it + 1;
   
     end
     
     
     % Update main variables
     hsave(:,is) = h;  
     Uxsave(:,is) = Ux;
     Uysave(:,is) = Uy;
     


surf(reshape(mesh.xv,Nx+1,Ny+1),reshape(mesh.yv,Nx+1,Ny+1),reshape(h,Nx+1,Ny+1));
drawnow
     

      
     
     
      
end




























