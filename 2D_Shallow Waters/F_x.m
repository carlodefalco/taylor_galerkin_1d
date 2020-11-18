function [Fhx, FU1x,FU2x] = F_x (h,Ux,Uy, g)

%   Fhx = zeros(size(h)); 
%   FU1x = zeros (size (h));
%   FU2x = zeros(size(h));
  
  % pure velocity component
  vx = Ux./h;
  vy = Uy./h;
  vy(h<=0) = 0;
  vx(h<=0) = 0;
  
  Fhx = Ux; 
  FU1x = Ux.*vx + g*(h.^2)/2; 
  FU2x = Ux.*vy;

end