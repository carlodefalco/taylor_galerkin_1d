function[dxh,dxUx,dxUy,dyh,dyUx,dyUy] = dU(mesh,h,Ux,Uy,dt)


dxh = zeros(mesh.NP,1);
dxUx = zeros(mesh.NP,1);
dxUy = zeros(mesh.NP,1);

dyh = zeros(mesh.NP,1);
dyUx = zeros(mesh.NP,1);
dyUy = zeros(mesh.NP,1);

for ip = 1:mesh.NP
    
    for j = 1:4
        for nq = 1:4
            
            dxh(ip) = dxh(ip) +...
               dt*h(mesh.conn(ip,j))*...
               (mesh.shgx(ip,nq,j))*mesh.wq(ip,nq);
           
           dyh(ip) = dyh(ip) +...
               dt*h(mesh.conn(ip,j))*...
               (mesh.shgy(ip,nq,j))*mesh.wq(ip,nq);
           
           
            dxUx(ip) = dxUx(ip) +...
               dt*Ux(mesh.conn(ip,j))*...
               (mesh.shgx(ip,nq,j))*mesh.wq(ip,nq);
           
            dyUx(ip) = dyUx(ip) +...
              dt*Ux(mesh.conn(ip,j))*...
               (mesh.shgy(ip,nq,j))*mesh.wq(ip,nq);
           
           
            dxUy(ip) = dxUy(ip) +...
              dt*Uy(mesh.conn(ip,j))*...
               (mesh.shgx(ip,nq,j)+mesh.shgy(ip,nq,j))*mesh.wq(ip,nq);
           
           dyUy(ip) = dyUy(ip) +...
               dt*Uy(mesh.conn(ip,j))*...
               (mesh.shgy(ip,nq,j))*mesh.wq(ip,nq);
           
        end
    end
end

return