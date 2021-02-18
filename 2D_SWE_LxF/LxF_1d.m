clear all;
close all;




tol = .001; 
x   = linspace (-400, 400, 50) .';
dx  = diff (x);

T   = 1000;
dt0 = .001;
dt  = .001;
tsave = linspace (0, T, 100);


h0 = 8 + sin (pi * x / 800);
u0 = zeros (size (x));

h = h0;
u = u0;

hsave = h;
usave = u;

t = 0;
tlist = t;
dtlist = dt;
errlist = 0;
for is = 2 : numel (tsave)
  while (t < tsave(is))
    tnew = t + dt;
    if (tnew > tsave(is))
      dt   = tsave(is) - t;
      tnew = tsave(is);
    end
    fprintf ('t = %g, dt = %g\n', tnew, dt)

    hnew = hstep (t, dt, x, dx, h, u);
    unew = ustep (t, dt, x, dx, h, u);

    hnews = hstep (t, dt/2, x, dx, h, u);
    unews = ustep (t, dt/2, x, dx, h, u);
    hnews = hstep (t+dt/2, dt/2, x, dx, hnews, unews);
    unews = ustep (t+dt/2, dt/2, x, dx, hnews, unews);
    

    herr = norm (hnew - hnews, inf) / dt;
    uerr = norm (unew - unews, inf) / dt;
    errmax = max (herr, uerr);
    
    if (errmax <= tol)
      t  = tnew;
      h  = hnew;
      u  = unew;
      dt = min (.5 * tol/errmax * dt, 1.5 * dt);
      tlist(end+1)  = t;
      dtlist(end+1) = dt;
      errlist(end+1)= errmax;
    else
      %keyboard
      dt = min (.5 * tol/errmax * dt, dt/2);
    end
    
    end

  hsave(:,is)  = h;
  usave(:, is) = u;
    end


for is = 1 : numel (tsave)

  figure (1)
  plot (x, hsave(:, is))
  axis ([-400, 400, 7, 9]);
  drawnow
  fprintf ("-dpng", "-S400,300", "-F:9", sprintf("h_%4.4d.png", is));
  
  figure (2)
  plot (x, usave(:, is))
  axis ([-400, 400, -10, 10]);
  drawnow
  fprintf ("-dpng", "-S400,300", "-F:9", sprintf("u_%4.4d.png", is));

end


% END OF THE CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xm = p0val (x)
  xm = .5 * (x(2:end) + x(1:end-1));
end

function Fh = fluxh (h, u)
  Fh = u;
end

function Sh = sourceh (h, u)
  Sh = 0 * h;
end

function Su = sourceu (h, u)
tau = .009;
  Su = - tau * u;
end

function Fu = fluxu (h, u)
 g = 9.81;
  Fu = u.^2 ./ h + (1/2) * g * h.^2;
end

function anew = astep (t, dt, x, dx, a, Flf, S)
V = 1.05;
  mass = .5 * ([dx; 0] + [0; dx]);
  
  % viscosita' artificiale (ovvero flusso correttivo) fissata
  da = V .* diff(a);
  
  % viscosita' artificiale (ovvero flusso correttivo) alla
  % Lax-Friedrichs
  % da = (dx / dt) .* diff(a);
  
  Fstar = Flf - da; 
  anew = a + dt * (([dx.*(Fstar .* (-1./dx)); 0] +...
                   [0; dx.*(Fstar .* (1./dx))])./ mass + S);
end

function hnew = hstep (t, dt, x, dx, h, u)
  Flf = p0val (fluxh (h, u));
  S   = sourceh (h, u);
  hnew = astep (t, dt, x, dx, h, Flf, S);
end

function unew = ustep (t, dt, x, dx, h, u)
  Flf = p0val (fluxu (h, u));
  S   = sourceu (h, u);
  unew = astep (t, dt, x, dx, u, Flf, S);
  unew([1 end]) = 0;
end

