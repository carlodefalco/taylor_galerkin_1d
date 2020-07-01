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
      
      u = L*x*exp(-t.^2);
      ux = L*exp(-t.^2);
      uxt = -2*L*t.*exp(-t.^2);
      ut = -2*L*t.*x.*exp(-t.^2);
      utt = -2*L*x.*exp(-t.^2)+4*t*L*x.*exp(-t.^2);
      s = ut + A*ux;
      st = utt + A*uxt;
    
  end
  
end
