function RHS = right_hand(msh, u, A, s,dt)

  d1 = zeros (msh.ndof, 1);
   d2 = zeros (msh.ndof, 1);
  for k = 1 : msh.nel
    for i = 1 : 2 
      for j = 1 : 2     
        d1(msh.conn(i, k)) = d1(msh.conn(i, k)) + ...
        - A * u(msh.conn(j, k)) * msh.shg(j, k) * msh.shp(i, i, k) * msh.h(k)/2 + ... %% advection
          s(msh.conn(i, k)) * msh.h(k) * msh.shp(i, j, k)  / 2;                      %% source
      
        d2(msh.conn(i, k)) = d2(msh.conn(i, k)) + ...
        - A^2 * u(msh.conn(j, k)) * msh.shg(j, k) * msh.shg(i, k) * msh.h(k);
      end
    end
    
  end
  
RHS = (dt)*d1 +1/2*(dt^2)*d2;


return


% This leads to:   s(i)*h/2 + A/h * h * (u(i+1)-u(i-1))/2