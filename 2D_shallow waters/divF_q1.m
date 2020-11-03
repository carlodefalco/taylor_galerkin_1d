function [dh_f,dUx_f,dUy_f] = divF_q1(mesh,Fhx_half,FU1x_half,FU2x_half,Fhy_half,FU1y_half,FU2y_half)


dh_f = zeros(mesh.NV,1);
dUx_f = zeros(mesh.NV,1);
dUy_f = zeros(mesh.NV,1);

for p = 1:mesh.NP
    
    dhloc = zeros(4,1);
    dUxloc = zeros(4,1);
    dUyloc = zeros(4,1);
    
    for itest = 1:4
        for jshape = 1:4
            for nq = 1:6
                
                dhloc(itest) = dhloc(itest)+...
                    (mesh.shg(1,nq,itest,p)*Fhx_half(p)+...
                    mesh.shg(2,nq,itest,p)*Fhy_half(p))*...
                    mesh.polygon(p).wq(nq);
                
                dUxloc(itest) = dUxloc(itest)+...
                    (mesh.shg(1,nq,itest,p)*FU1x_half(p)+...
                    mesh.shg(2,nq,itest,p)*FU1y_half(p))*...
                    mesh.polygon(p).wq(nq);
                
                dUyloc(itest) = dUyloc(itest)+...
                    (mesh.shg(1,nq,itest,p)*FU2x_half(p)+...
                    mesh.shg(2,nq,itest,p)*FU2y_half(p))*...
                    mesh.polygon(p).wq(nq);
                
            end
        end
        
    end
            dh_f(mesh.conn(p,1:4)) = dh_f(mesh.conn(p,1:4)) + dhloc;
            dUx_f(mesh.conn(p,1:4)) = dUx_f(mesh.conn(p,1:4)) + dUxloc; 
            dUy_f(mesh.conn(p,1:4)) = dUy_f(mesh.conn(p,1:4)) + dUyloc;     
       
        
    
    
end


                
                
                