function [Sh, SU1,SU2] = source (h, Z_x, Z_y, g)

  Sh = zeros(size(h));
  SU1 = zeros (size (h));
  SU2 = zeros (size (h));
%   v = U ./ h;
%   v(h <= 0) = 0;
  
  Sh(:)   = 0;
  SU1      = -g*h.*Z_x;
  SU2      = -g*h.*Z_y;
  
end