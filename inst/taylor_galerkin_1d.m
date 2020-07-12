clear all
close all

coefficients

msh.ndof       = N+1;
msh.nel        = msh.ndof - 1;
msh.x          = linspace (0, L, msh.ndof).';
msh.conn       = [1:msh.ndof-1; 2:msh.ndof];
msh.h          = diff (msh.x).';
msh.shg(1, :)  = -1./msh.h;
msh.shg(2, :)  = +1./msh.h;
msh.shp(1,1,:) = ones (msh.nel, 1);
msh.shp(2,2,:) = ones (msh.nel, 1);
t              = linspace (0, T, M+1);

mass      = h * ones (size (msh.x));
mass(1)   =  mass(1) / 2;
mass(end) =  mass(end) / 2;

source   = @example_source;
solution = @example_solution;
u(:, 1)  = solution (0, msh.x);
w(:,1) = solution(0,msh.x);

for n = 1 : M

  [uex, ~, uext] = solution (t(n+1), msh.x);
  U1   = uex(1);
  UN1  = uex(end);
  [s,st,sx]  = source (t(n), msh.x);
  
  S = mass.*s;
  St = mass.*st;
  Sx = mass.*sx;

% S = source_2(msh, s);
% St = source_2(msh, st);
% Sx = source_2(msh, sx);
  
  

  u(1, n+1)   = U1;
  w(1,n+1) = U1;
  
  du1 = tg2d1 (msh, u(:,n), A, s);
  du2 = tg2d2 (msh, u(:,n), A, s);
  
  
  du1 = du1 ./ mass;
  du2 = du2 ./ mass;
  
  

  w(2:end,n) = u(2:end,n) + 1/3*dt*(s(2:end)+du1(2:end)) + 1/9*(dt^2)*(St(2:end)+A*Sx(2:end)+du2(2:end));
  
  dw2 = tg2d2 (msh, w(:,n), A, s);
  dw2 = dw2./mass;
  
  u(2:end,n+1) = u(2:end,n) + dt*(s(2:end)+du1(2:end)) + (dt^2/2)*(St(2:end)+A*Sx(2:end)+dw2(2:end));
  
  
  
  
  
  
  
  
  
  err(n+1) = trapz (msh.x, abs (uex- u(:, n+1)).^2);
  
  
  figure (1)
  subplot (3, 1, 1)
  plot (msh.x, uex,'b--', msh.x, u(:, n+1), 'r--','linewidth',1.2);
  legend ('exact', 'computed')
  
  subplot (3, 1, 2)
  plot (msh.x, abs (uex- u(:, n+1)));
  title ('error in space')

  subplot (3, 1, 3)
  plot (t(1:n+1), err);
  title ('L^2-norm of error in space vs time')
  drawnow
  
end


disp(max(err));