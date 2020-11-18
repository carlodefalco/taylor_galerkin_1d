function [dh_f,dUx_f,dUy_f] = divF_q1(mesh,Fhx_half,FU1x_half,FU2x_half,Fhy_half,FU1y_half,FU2y_half)


dh_f = zeros(mesh.NV,1);
dUx_f = zeros(mesh.NV,1);
dUy_f = zeros(mesh.NV,1);

for p = 1:mesh.NP
    

 Qe = mesh.conn(p,:);
    
    dhloc = zeros(4,1);
    dUxloc = zeros(4,1);
    dUyloc = zeros(4,1);
    
    for itest = 1:4
        for jshape = 1:4
            for nq = 1:4
                
                dhloc(itest) = dhloc(itest)+...
                    (mesh.shgx(p,nq,itest)*Fhx_half(p)+...
                    mesh.shgy(p,nq,itest)*Fhy_half(p))*...
                    mesh.wq(p,nq);
                
                dUxloc(itest) = dUxloc(itest)+...
                    (mesh.shgx(p,nq,itest)*FU1x_half(p)+...
                    mesh.shgy(p,nq,itest)*FU1y_half(p))*...
                    mesh.wq(p,nq);
                
                dUyloc(itest) = dUyloc(itest)+...
                  (mesh.shgx(p,nq,itest)*FU2x_half(p)+...
                    mesh.shgy(p,nq,itest)*FU2y_half(p))*...
                    mesh.wq(p,nq);
                
            end
        end
        
    end
            dh_f(Qe) = dh_f(Qe) + dhloc;
            dUx_f(Qe) = dUx_f(Qe) + dUxloc; 
            dUy_f(Qe) = dUy_f(Qe) + dUyloc;     
       
        
    
    
end


                
                
                