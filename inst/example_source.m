function [s, st] = example_source (t, x)

  coefficients

  if (testnum == 1)
    xi  = t/T;
    xit = 1/T;
    u   =  ((x-L).^2/(L^2/(1-xi)) + xi);
    ux  =  (2*(1 - t/T) .* (-L + x)) / L^2;
    uxt =  (2*(- 1/T) .* (-L + x)) / L^2;
    ut  =  1/T - (-L + x).^2 ./(L.^2 * T);
    utt =  0;
    s   = ut + A* ux;
    st  = utt + A* uxt;
  elseif (testnum == 2)
    s = zeros (size (x));
    st = zeros (size (x));
    
   elseif (testnum == 3)
    u = t/T.*exp(sin(L*x));
    ut = 1/T.*exp(sin(L*x));
    ux = t/T*L*cos(L*x).*exp(sin(L*x));
    uxt = 1/T.*L.*cos(L*x).*exp(sin(L*x));
    utt = 0;
    s = ut + A*ux;
    st = utt + A*uxt;
    
  elseif (testnum == 4)
    u = 3*sin(x.^2);
    ut = zeros(size(x));
    ux = 6*x.*cos(x.^2);
    uxt = zeros(size(x));
    utt = zeros(size(x));
    s = ut + A*ux;
    st = utt + A*uxt;
    sx = 6*cos(x.^2) - 12*x.^2.*sin(x.^2); 
    
    
    
    
    
    
    
  end
  
end
