clear all
close all

    
coefficients

msh.ndof       = N+1;
msh.nel        = msh.ndof - 1;
msh.x          = linspace(0, L, msh.ndof).';


msh.conn       = [1:msh.ndof-1; 2:msh.ndof];
msh.h          = diff (msh.x).';
msh.shg(1, :)  = -1./msh.h;
msh.shg(2, :)  = +1./msh.h;
msh.shp(1,1,:) = ones (msh.nel, 1);
msh.shp(2,2,:) = ones (msh.nel, 1);
t              = linspace (0, T, M+1);



mass      = h * ones (size (msh.x));

mass(1)   =  mass(1) / 2;
mass(end) =  mass(end) / 2;

source   = @example_source;
solution = @example_solution;
u(:, 1)  = solution (0, msh.x);






for n = 1 : M

  [uex, ~, uext] = solution (t(n+1), msh.x);
  U1   = uex(1);
  UN1  = uex(end);

  
  u(1, n+1)   = U1;
  u(n+1, n+1) = UN1;
  
  s = source (t(n), msh.x); % it was n+1
  s1 = source_2(msh, s);
  
  
  du1 = tg3d1 (msh, u(:,n), A, s);
  du2 = tg3d2 (msh, u(:,n), A, s);
  lhs = tg3lhs (msh, u(:,n), A, s, dt);
  
  [s,st,sx]  = source (t(n), msh.x);
  stN = source_2(msh, st);
  sxN = source_2(msh,sx);
  
  
  [s,st,sx]  = source (t(n+1), msh.x);
  stN1 = source_2(msh, st);
  sxN1 = source_2(msh,sx);
  

   
  
  rhs([1 msh.ndof],1) = [U1 UN1];




  rhs(2:end-1,1) = mass(2:end-1).*u(2:end-1, n) + (dt)*(du1(2:end-1)+ s1(2:end-1)) + ...
                   (dt^2/2)*du2(2:end-1) + (dt^2/3)*stN(2:end-1) + (dt^2/6)*stN1(2:end-1) + ...
                   (dt^2/3)*(A*sxN(2:end-1)) + (dt^2)/6*(A*sxN1(2:end-1));


  
  bignumber = 1e13 ;
  lhs(1,1) =  bignumber;
  lhs(end,end) =  bignumber;
  rhs(1) = U1*bignumber;
  rhs(end) = UN1*bignumber;
  u(:, n+1)   = lhs \ rhs;



  err(n+1) = trapz (msh.x, abs (uex- u(:, n+1)).^2);
  
  
  
  subplot (3, 1, 1)
  plot (msh.x, uex, 'b--', msh.x, u(:, n+1), 'r--');
  legend ('exact', 'computed')
  
                                
  
  subplot (3, 1, 2)
  plot (msh.x, abs(uex- u(:, n+1)));
  title ('error in space')

                           
  
  subplot (3, 1, 3)
  plot (t(1:n+1), err);
  title ('L^2-norm of error in space vs time')
  drawnow

end

                              


