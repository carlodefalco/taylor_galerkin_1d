function [s_N] = source_2(msh, source)


  s_N = zeros (msh.ndof, 1);
  for k = 1 : msh.nel
    for i = 1 : 2
      for j = 1 : 2
        s_N(msh.conn(i, k)) = s_N(msh.conn(i, k)) + ...
        source(msh.conn(i, k)) * msh.h(k) * msh.shp(i, j, k)/2;
    
      end
    end
  end
  
end