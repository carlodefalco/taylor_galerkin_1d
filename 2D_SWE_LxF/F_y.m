function [Fhy, FU1y,FU2y] = F_y (h, Ux,Uy, g)

%   Fhy = zeros(size(h)); 
%   FU1y = zeros (size (h));
%   FU2y = zeros(size(h));
  
  % pure velocity component
  vx = Ux./h;
  vy = Uy./h;
  vy(h<=0) = 0;
  vx(h<=0) = 0;
  
  Fhy = Uy; 
  FU1y = Uy.*vx; 
  FU2y = Uy .* vy + g * (h .^2) / 2;

end