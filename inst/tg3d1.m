function d1 = tg3d1 (msh, u, A, s)

  d1 = zeros (msh.ndof, 1);
  for k = 1 : msh.nel
    for i = 1 : 2
      for j = 1 : 2     
        d1(msh.conn(i, k)) = d1(msh.conn(i, k)) + ...
        - A * u(msh.conn(j, k)) * msh.shg(j, k) * msh.shp(i, i, k) * msh.h(k)/2 + ... %% advection
          
      end
    end
  end
  
end

