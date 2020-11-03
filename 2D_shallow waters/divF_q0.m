function [dh,dUx,dUy] = divF_q0(mesh,Fhx,FU1x,FU2x,Fhy,FU1y,FU2y)

dh = zeros(mesh.NP,1);
dUx = zeros(mesh.NP,1);
dUy = zeros(mesh.NP,1);

for p = 1:mesh.NP
    
    hloc = zeros(4,p);
    Uxloc = zeros(4,p);
    Uyloc = zeros(4,p);
    
    for itest = 1:4
        for jshape = 1:4
            for nq = 1:6
                
                hloc(itest,p) = hloc(itest,p)+...
                    -mesh.area(p)*...
                    (mesh.shg(1,nq,jshape,p)*Fhx(mesh.conn(p,jshape))+...
                    mesh.shg(2,nq,jshape,p)*Fhy(mesh.conn(p,jshape)))*...
                    mesh.polygon(p).wq(nq);
                
                

                Uxloc(itest,p) = Uxloc(itest,p)+...
                    -mesh.area(p)*...
                    (mesh.shg(1,nq,jshape,p)*FU1x(mesh.conn(p,jshape))+...
                    mesh.shg(2,nq,jshape,p)*FU1y(mesh.conn(p,jshape)))*...
                    mesh.polygon(p).wq(nq);
                
                Uyloc(itest,p) = Uyloc(itest,p)+...
                    -mesh.area(p)*...
                    (mesh.shg(1,nq,jshape,p)*FU2x(mesh.conn(p,jshape))+...
                    mesh.shg(2,nq,jshape,p)*FU2y(mesh.conn(p,jshape)))*...
                    mesh.polygon(p).wq(nq);
                
            end
        end
        
    end
        
        dh(p) = sum(hloc(1:4,p))/4; 
        dUx(p) = sum(Uxloc(1:4,p))/4;
        dUy(p) = sum(Uyloc(1:4,p))/4;
        
    
        
end

                    
                
                
                
                
                
                
                
                
                