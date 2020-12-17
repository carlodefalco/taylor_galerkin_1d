function L = Mass(mesh)

% This script return the mass matrix using the Trapezium Quadrature rule.
% In this way we have a diagonal matrix (lumped).

L = sparse(mesh.NV,mesh.NV);

for ip = 1:mesh.NP
    
    Qe = mesh.conn(ip,:);
    Mloc = zeros(4);
    
for i = 1:4
    for j = 1:4
        for nq = 1:4
            
            Mloc(i,j) = Mloc(i,j) +...
                mesh.shape(ip,nq,i)*mesh.shape(ip,nq,j)*mesh.wq(ip,nq);
        end
    end
end

    
    
    L(Qe,Qe) = L(Qe,Qe) + Mloc;
    
    clear Mloc;
    
end

return