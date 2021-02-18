function mesh = meshinfo(mesh)

%{
 This function return some infos about the mesh used.
All the informations are located in the struct MESH and these are:

* MESH.CONN --> Connectivity matrix
* MESH.X and MESH.Y --> coordinates of all dofs
* MESH.AREA --> area of each element
* MESH.SHP --> shape functions evaluated on each element
* MESH.SHG --> shape functions gradient
%}


% Load Connectivity Matrix
for p = 1:mesh.NP
    mesh.conn(p,1:4) = mesh.polygon(p).vertices(:);  
end

% Load mesh dofs (vertices)
for v = 1:mesh.NV
mesh.x(v) = mesh.vertex(v).x;
mesh.y(v) = mesh.vertex(v).y;
end




for ip = 1:mesh.NP
    
    % Storing the surface of each polygon
     mesh.area(ip,1) = mesh.polygon(ip).area;
    
 x1 =   mesh.vertex(mesh.polygon(ip).vertices(1)).x ;
 x2 =   mesh.vertex(mesh.polygon(ip).vertices(2)).x ;
 x3 =   mesh.vertex(mesh.polygon(ip).vertices(3)).x ;
 x4 =   mesh.vertex(mesh.polygon(ip).vertices(4)).x ;
 
 y1 =   mesh.vertex(mesh.polygon(ip).vertices(1)).y ;
 y2 =   mesh.vertex(mesh.polygon(ip).vertices(2)).y ;
 y3 =   mesh.vertex(mesh.polygon(ip).vertices(3)).y ;
 y4 =   mesh.vertex(mesh.polygon(ip).vertices(4)).y ;
 
 mesh.xv(ip,1:4) = [x1 x2 x3 x4];
 mesh.yv(ip,1:4) = [y1 y2 y3 y4];
 mesh.wq(ip,1:4) = .25*[mesh.area(ip) mesh.area(ip) mesh.area(ip) mesh.area(ip)];
 

    
    mesh.hx(ip,1) = mesh.xv(ip,2)-mesh.xv(ip,1);
    mesh.hy(ip,1) = mesh.yv(ip,3)-mesh.yv(ip,1);
    
    mesh.shape(ip,1:4,1) = [1 0 0 0 ];
    mesh.shape(ip,1:4,2) = [0 1 0 0 ];
    mesh.shape(ip,1:4,3) = [0 0 1 0 ];
    mesh.shape(ip,1:4,4) = [0 0 0 1 ];
    
    mesh.shgx(ip,1:4,1) = [-1/mesh.hx(ip,1) -1/mesh.hx(ip,1) 0 0];
    mesh.shgx(ip,1:4,2) = [1/mesh.hx(ip,1) 1/mesh.hx(ip,1) 0 0];
    mesh.shgx(ip,1:4,3) = [0  0  1/mesh.hx(ip,1) 1/mesh.hx(ip,1)];
    mesh.shgx(ip,1:4,4) = [0  0  -1/mesh.hx(ip,1)  -1/mesh.hx(ip,1)];
    
    mesh.shgy(ip,1:4,1) = [-1/mesh.hy(ip,1) 0  0 -1/mesh.hy(ip,1)  ];
    mesh.shgy(ip,1:4,2) = [ 0 -1/mesh.hy(ip,1) -1/mesh.hy(ip,1) 0];
    mesh.shgy(ip,1:4,3) = [ 0 1/mesh.hy(ip,1) 1/mesh.hy(ip,1) 0];
    mesh.shgy(ip,1:4,4) = [  1/mesh.hy(ip,1) 0  0 1/mesh.hy(ip,1) ];


    
    
    
    
end

   
      
      




