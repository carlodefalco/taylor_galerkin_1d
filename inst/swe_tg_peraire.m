1;

function x = jacobi (A, P, b, nit = 5, damp = 1)
  
  res = b; 
  x   = damp*(P \ res);
  for it = 2 : nit
    x = x + damp*(P \ res);
    res = (b - A * x);
  endfor
endfunction

function [Fh, FU] = flux (h, U, g = 9.81)

  Fh = FU = zeros (size (h));
  v = U ./ h;
  v(h <= 0) = 0;
  
  Fh = U;
  FU = U .* v + g * (h .^2) / 2;   

endfunction


function [Sh, SU] = source (h, U, gradz, g = 9.81)

  Sh = SU = zeros (size (h));
  v = U ./ h;
  v(h <= 0) = 0;
  
  Sh(:)   = 0;
  SU      = -g * h .* gradz;
  
endfunction

L              = 400;
N              = 50;
msh.ndof       = N+1;
msh.nel        = msh.ndof - 1;
msh.x          = linspace (-L, L, msh.ndof).';
msh.xm         = (msh.x(2:end)+msh.x(1:end-1))/2;
msh.conn       = [1:msh.ndof-1; 2:msh.ndof];
msh.h          = diff (msh.x);
msh.shg(1, :)  = -1./msh.h;
msh.shg(2, :)  = +1./msh.h;


mass           =  sparse (diag (([msh.h;0] + [0;msh.h])/3) +
                          diag (msh.h/6, 1) + diag (msh.h/6, -1));
Dmass          = diag (sum (mass, 2));
%%Dmass = diag (diag (mass));

gradz(1:msh.ndof,1) = 0;%-2/(2*L);
g = 9.81;
Z = cumtrapz (msh.x, gradz);
Zm = (Z(2:end)+Z(1:end-1))/2;
h = h0 = 8 - sin (pi*msh.x/2/L) - Z;% 5 - Z;
U = U0 = zeros (size (msh.x));

T    = 1500;
t    = 0;
it   = 1;
dt   = 1e-2;
Njac = 100;

nsave = round (T/5);
tsave = linspace (0, T, nsave);

hsave = h0;
Usave = U0;

figure (1)
plot (msh.x, h , 'linewidth', 1,
      msh.x, Z , 'linewidth', .5)
grid on
title (sprintf ("t = %g", tsave(1)))
drawnow
figure (2)
plot (msh.x, U, 'k-', 'linewidth', .5)
title (sprintf ("t = %g", tsave(1)))
drawnow
  
for is = 2 : nsave 
  while (t(it) < tsave (is))

    t(it+1) = t(it) + dt;
    if (t(it+1) > tsave (is))
      t(it+1) = tsave (is);
      dt = t(it+1) - t(it);
    endif
    printf ("it = %d, t = %g, dt = %g\n", it, t(it+1), dt)

   
    [Fh, FU] = flux (h, U, g);    
    hm = (h(2:end)+h(1:end-1))/2;
    Um = (U(2:end)+U(1:end-1))/2;
    [~,  SU] = source (h, U, gradz, g);

    dh = zeros (msh.nel, 1);
    dU = zeros (msh.nel, 1);
    for k = 1 : msh.nel
      dh(k) = dh(k) + ...
      -(Fh(msh.conn(2, k)) - Fh(msh.conn(1, k)));
    
      dU(k) = dU(k) + ...
      -(FU(msh.conn(2, k)) - FU(msh.conn(1, k))) + ...
      (SU(msh.conn(2, k)) + SU(msh.conn(1, k))) * msh.h(k)/2;
      
    end 
    
    wh = hm + (dt/2) * (dh./msh.h(:));
    wU = Um + (dt/2) * (dU./msh.h(:));

    dh = zeros (msh.ndof, 1);
    dU = zeros (msh.ndof, 1);
    [Fh, FU] = flux (wh, wU, g);
    [~,  SU] = source (wh, wU, (gradz(1:end-1)+gradz(2:end))/2, g);

    for k = 1 : msh.nel
      for i = 1 : 2
        dh(msh.conn(i, k)) = dh(msh.conn(i, k)) + ...
                             Fh(k) * msh.shg(i, k) * msh.h(k);

         dU(msh.conn(i, k)) = dU(msh.conn(i, k)) + ...
                              FU(k) * msh.shg(i, k) * msh.h(k) + ...
                              SU(k) * msh.h(k)/2;
      end
    end
    
    h = h + dt * (jacobi (mass, Dmass, dh, Njac));
    U(2:end-1) = U(2:end-1) + dt * (jacobi (mass(2:end-1,2:end-1), Dmass(2:end-1,2:end-1), dU(2:end-1), Njac));
    U([1, end]) = 0;
    v = U ./ h;
    v(h <= 0) = 0;
    
    dtold = dt;
    dt = (.09/sqrt(3)) * min (min (msh.h) ./ (norm (v, inf) + norm (sqrt (g * h), inf)));
    dt = min (dt, 2 * dtold);        
    
    ++it;
  endwhile
 
  figure (1)
  plot (msh.x, h+Z, 'linewidth', 1,
        msh.xm, wh+Zm , 'linewidth', 1,
        msh.x, Z , 'linewidth', .5)
  grid on
  axis ([-L, L, 6, 10])
  title (sprintf ("t = %g", tsave(is)))
  drawnow
  figure (2)
  plot (msh.x, U, 'k-', 'linewidth', .5,
        msh.xm, wU , 'linewidth', .5)
  title (sprintf ("t = %g", tsave(is)))
  drawnow

  hsave (:, is) = h;
  Usave (:, is) = U;

  endfor
