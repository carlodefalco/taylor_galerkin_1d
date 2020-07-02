function d1 = tg3d2 (msh, u, A, s)

  d1 = zeros (msh.ndof, 1);
  for k = 1 : msh.nel
    for i = 1 : 2
      for j = 1 : 2
        d1(msh.conn(i, k)) = d1(msh.conn(i, k)) + ...
        - (2/3) * A^2 * u(msh.conn(j, k)) * msh.shg(j, k) * msh.shg(i, k) * msh.h(k);
      end
    end
  end
  
end


