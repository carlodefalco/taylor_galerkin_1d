function [s, st] = example_source (t, x)

  coefficients
  xi  = t/T;
  xit = 1/T;
  u   =  ((x-L).^2/(L^2/(1-xi)) + xi);
  ux  =  (2*(1 - t/T) .* (-L + x)) / L^2;
  uxt =  (2*(- 1/T) .* (-L + x)) / L^2;
  ut  =  1/T - (-L + x).^2 ./(L.^2 * T);
  utt =  0;
  s   = ut + A* ux;
  st  = utt + A* uxt;
  
endfunction
