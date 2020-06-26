clear all
%close all


coefficients
x = linspace (0, L, N+1)(:);
t = linspace (0, T, M+1)(:);

mass = h * ones (size (x));
mass(1) =  mass(1) / 2;
mass(end) =  mass(end) / 2;

source  = @example_source;
u(:, 1) = example_solution (0, x);



for n = 1 : M

  [uex, ~, uext] = example_solution (t(n+1), x);
  U1   = uex(1);
  UN1  = uex(end);
  s  = source (t(n+1), x);
  
  du1  = zeros (size (x));
  du2  = zeros (size (x));

  du1([1,N+1])+=   s([1,N+1]) * h/2;
  du1(2:N)    +=   s(2:N) * h;
  
  du1(1:N)    += - A * u(2:end,n)   / 2;
  du1(2:N+1)  +=   A * u(1:end-1,n) / 2;

  
  du2(1:N+1)  += - 2* A^2 * u(:,n)    / h;
  du2(2:N+1)  +=   A^2 * u(1:end-1,n) / h;
  du2(1:N)    +=   A^2 * u(2:end,n)   / h;

  du1 = du1 ./ mass;
  du2 = du2 ./ mass;
  
  
  u(:, n+1) = u(:, n) + (dt) * du1 + (dt)^2/2 * du2;
  u(1, n+1) = U1;
  u(end, n+1) = UN1;
  
  figure (1)
  subplot (3, 1, 1)
  plot (x, uex, x, u(:, n+1), 'x-');
  legend ('exact', 'computed')
  drawnow

  figure (1)
  subplot (3, 1, 2)
  plot (x, abs (uex- u(:, n+1)));
  title "error in space"
  drawnow

  figure (1)
  subplot (3, 1, 3)
  err(n+1) = trapz (x, abs (uex- u(:, n+1)).^2);
  plot (t(1:n+1), err);
  title "2-norm of error in space vs time"
  drawnow
endfor 

