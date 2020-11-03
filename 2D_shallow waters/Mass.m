function M = Mass(mesh)

M = sparse(mesh.NV,mesh.NV);

for ip = 1:mesh.NP
    
    Mloc = zeros(4);
    
    for in = 1:4
        for jn = 1:4
            
            Mloc(in,jn) = Mloc(in,jn) + ...
            sum(mesh.shp(:,in,ip).*mesh.shp(:,jn,ip).*(mesh.polygon(ip).wq(:)));
            
        end
        
    end
    
    M(mesh.conn(ip,1:4),mesh.conn(ip,1:4)) = M(mesh.conn(ip,1:4),mesh.conn(ip,1:4)) + Mloc;
    
end

return

        
        