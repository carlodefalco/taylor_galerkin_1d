function [st_N] = source_2(msh, s)


  st_N = zeros (msh.ndof, 1);
  for k = 1 : msh.nel
    for i = 1 : 2
      for j = 1 : 2
        st_N(msh.conn(i, k)) = st_N(msh.conn(i, k)) + ...
        s(msh.conn(i, k)) * msh.h(k) * msh.shp(i, j, k)/2;
    
      end
    end
  end
  
end