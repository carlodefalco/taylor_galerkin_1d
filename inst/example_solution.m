function [u, ux, ut] = example_solution (t, x)

  coefficients
  
  if (testnum == 1)
    xi  = t/T;
    xit = 1/T;
    u   =  ((x-L).^2/(L^2/(1-xi)) + xi);
    ux  =  (2*(1 - t/T) .* (-L + x)) / L^2;
    ut  =  1/T - (-L + x).^2 ./(L.^2 * T);
  elseif (testnum == 2)
    kappa = (omega/A)*(L/T);
    u   = sin ( 2*pi * (kappa*x/L-omega*t/T) );
    ut  = -2*pi*omega/T * cos ( 2*pi * (kappa*x/L-omega*t/T) );
    ux  = 2*pi*kappa/L * cos ( 2*pi * (kappa*x/L-omega*t/T) );
  endif
  
endfunction
