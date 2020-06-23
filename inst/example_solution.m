function [u, ux, ut] = example_solution (t, x)

  coefficients
  xi  = t/T;
  xit = 1/T;
  u   =  ((x-L).^2/(L^2/(1-xi)) + xi);
  ux  =  (2*(1 - t/T) .* (-L + x)) / L^2;
  ut  =  1/T - (-L + x).^2 ./(L.^2 * T);
endfunction
