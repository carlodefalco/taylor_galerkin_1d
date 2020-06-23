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
  
  [s, st] = source (t(n), x);    
  du1  = zeros (size (x));
  du2  = zeros (size (x));

  du1(1:N)    += s(1:N) * h/2;
  du1(2:N+1)  += s(2:N+1) * h/2;
  ux           = (u(2:end,n) - u(1:end-1,n))/h;
  wxl          = -1/h;
  wxr          =  1/h;
  du1(1:N)    += - A * ux * h/2;
  du1(2:N+1)  += - A * ux * h/2;
  du2(1:N)    += - A^2 * ux * wxl * h/2;
  du2(1:N)    +=   A * s(1:N) * wxl * h/2;
  du2(1:N)    +=   st(1:N) * h/2;
  du2(2:N+1)  += - A^2 * ux * wxr * h/2;
  du2(2:N+1)  +=   A * s(2:N+1) * wxr * h/2;
  du2(2:N+1)  +=   st(2:N+1) * h/2;
  
  du1 = du1 ./ mass;
  du1(1) = 0;
    
  du2 = du2 ./ mass;
  du2(1) = 0;   

  dd = j / 1e2;
  
  u(:, n+1) = max (0, u(:, n) + dt * du1 + dt^2/2 * du2);
  u(1, n+1) = U1;

  figure (1)
  plot (x, uex, x, u(:, n+1), 'x-');
  drawnow

  figure (2)
  plot (x, abs (uex- u(:, n+1)));
  drawnow

  figure (3)
  err(n+1) = trapz (x, abs (uex- u(:, n+1)).^2);
  plot (t(1:n+1), err);
endfor 

