function lhs = tg3lhs (msh, u, A, s, dt)

  lhs = sparse (msh.ndof, msh.ndof);
  for k = 1 : msh.nel
    for i = 1 : 2
      lhs(msh.conn(i, k), msh.conn(i, k)) = lhs(msh.conn(i, k), msh.conn(1, k)) + msh.h(k) / 2;                                            
      for j = 1 : 2
        lhs(msh.conn(i, k), msh.conn(j, k)) = lhs(msh.conn(i, k), msh.conn(j, k))+...
        + (dt^2*A^2)/6 * msh.shg(j, k) * msh.shg(i, k) * msh.h(k); 
      end
    end
  end
  
end


