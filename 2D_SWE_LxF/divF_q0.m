function [dh,dUx,dUy] = divF_q0(mesh,Fhx,FU1x,FU2x,Fhy,FU1y,FU2y)

dh = zeros(mesh.NP,1);
dUx = zeros(mesh.NP,1);
dUy = zeros(mesh.NP,1);

      for p = 1:mesh.NP
   
    
    
        for j = 1:4
            for nq = 1:4
                
                dh(p) = dh(p)+...
                    -(mesh.shgx(p,nq,j)*Fhx(mesh.conn(p,j))+...
                    mesh.shgy(p,nq,j)*Fhy(mesh.conn(p,j)))*...
                    mesh.wq(p,nq);
                
                

                dUx(p) = dUx(p)+...
                    -(mesh.shgx(p,nq,j)*FU1x(mesh.conn(p,j))+...
                    mesh.shgy(p,nq,j)*FU1y(mesh.conn(p,j)))*...
                    mesh.wq(p,nq);
                
                dUy(p) = dUy(p)+...
                    -(mesh.shgx(p,nq,j)*FU2x(mesh.conn(p,j))+...
                    mesh.shgy(p,nq,j)*FU2y(mesh.conn(p,j)))*...
                    mesh.wq(p,nq);
                
            end
        end
        
     end
        

        
    
        
return
                    
                
                
                
                
                
                
                
                
                