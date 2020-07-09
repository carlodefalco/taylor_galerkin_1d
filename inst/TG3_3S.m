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
wt(:,1) = solution(0,msh.x);
wtt(:,1) = solution(0,msh.x);






for n = 1 : M

  [uex, ~, uext] = solution (t(n+1), msh.x);
  U1   = uex(1);
  UN1  = uex(end);

  
  u(1, n+1)   = U1;
  u(n+1, n+1) = UN1;
  
  wt(1,n+1) = U1;
 wt(n+1,n+1) = UN1;
  
  wtt(1,n+1) = U1;
 wtt(n+1,n+1) = UN1;
  
 [ s,~,~,sxx] = source (t(n), msh.x); % 
  
     s = mass.*s;
     sxx = mass.*sxx;



  du1 = tg3d1 (msh, u(:,n), A, s);
  du2 = tg3d2 (msh, u(:,n), A, s);
  
  du1 = du1./mass;
  du2 = du2./mass;
  

 % 1-STEP:    wt^n = M*u^n - dt/3*(a*Delta(u^n)) - dt/3*M*s 
 % where:
 %              M = (1+1/9*d^2_x) <---- d^2 is approximated by tg3d2 (approx.
 %              of the II derivative)
 %
 %              Delta(u) <--------- is approximated by tg3d1 (approx. of
 %              the I derivative)
 %
 %              M*s = s + M*s =  s + 1/9*sxx <-------  s,sxx are known functions
 % 
 

 wt(2:end,n) = u(2:end,n) + 1/9*du2(2:end) - (dt/3)*A*du1(2:end) -dt/3*(s(2:end)+ sxx(2:end)) ;
 
 % 2-STEP:    wtt^n = M*u^n - dt/2*(a*Delta(wt^n)) - dt/2*M*s
 % 
 %            wt1 = Delta(u^n_tilde) is approximated by tg3d1 evaluated in wt
 %            instead of u^n
 
 wt1 = tg3d1(msh,wt(:,n),A,s);
 
 wt1 = wt1./mass;
 
 wtt(2:end,n) = u(2:end,n) + 1/9*du2(2:end) - (dt/2)*A*wt1(2:end) -dt/2*(s(2:end)+ sxx(2:end)) ;
 
 
 % 3-STEP:    u^n+1 = M*u^n - dt*(a*Delta(wtt^n)) - dt/2*M*s
 
 
 wtt1 = tg3d1(msh,wtt(:,n),A,s);
 
 wtt1 = wtt1./mass;
 
 u(2:end,n+1) = u(2:end,n) + 1/9*du2(2:end) - dt*A*wtt1(2:end) -dt*(s(2:end)+ sxx(2:end))  ; 
 
 
 
  
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

                              
disp(max(err))

