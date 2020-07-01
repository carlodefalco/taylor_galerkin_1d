clear all
close all

coefficients

msh.ndof       = N+1; % n of dofs
msh.nel        = msh.ndof - 1; % n of elements
msh.x          = linspace (0, L, msh.ndof).';
msh.conn       = [1:msh.ndof-1; 2:msh.ndof]; % Connectivity matrix
msh.h          = diff (msh.x).'; % space step
msh.shg(1, :)  = -1./msh.h;  
msh.shg(2, :)  = +1./msh.h;
msh.shp(1,1,:) = ones (msh.nel, 1);  % shape fun phi_k
msh.shp(2,2,:) = ones (msh.nel, 1);
t              = linspace(0, T, M+1);

% mass matrix
mass = h*eye(size(msh.x,1));
mass(1,1)   =  mass(1,1) / 2;
mass(end,end) =  mass(end,end) / 2;

% constant h, linear fem Stiffness matrix - Only for one dimension 
Kh = 1/h*[diag(2*ones(1,size(msh.x,1)))+diag(-1*ones(1,size(msh.x,1)-1),1)+...
    diag(-1*ones(1,size(msh.x,1)-1),-1)];
Kh(end,end) = 1/h;
Kh = A^2*(dt^2/6)*Kh;

% LHS of the linear system D = M + a^2*(dt)^2/6*Kh
D = mass+Kh; 

source   = @example_source;
solution = @example_solution;
u(:, 1)  = solution (0, msh.x);



for n = 1 : M

  [uex, ~, uext] = solution (t(n+1), msh.x);
  U1   = uex(1);
  UN1  = uex(end);
  s  = source (t(n+1), msh.x);

  u(1, n+1)   = U1;
  
  
  
 RHS = right_hand(msh, u(:,n), A, s,dt);
 
 RHS = D\RHS;
 

 
u(2:end,n+1) = u(2:end,n) + RHS(2:end);


  
  
  

  err(n+1) = trapz (msh.x, abs (uex- u(:, n+1)).^2);
  figure (1)
  subplot (3, 1, 1)
  plot (msh.x, uex,'b--', msh.x, u(:, n+1), 'r--');
  legend ('exact', 'computed')
  
  subplot (3, 1, 2)
  plot (msh.x, abs (uex- u(:, n+1)));
  title ('error in space');

  subplot (3, 1, 3)
  plot (t(1:n+1), err);
  title ('L^2-norm of error in space vs time')
  drawnow
  
end


