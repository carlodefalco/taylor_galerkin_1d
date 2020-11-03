
function mesh = meshinfo(mesh)

%{
 This function return some infos about the mesh used.
All the informations are located in the struct MESH and these are:

* MESH.CONN --> Connectivity matrix
* MESH.XV and MESH.YV --> coordinates of all dofs (only vertices for Q1 elements!)
* MESH.POLYGON(:).XQ/YQ/WQ --> Vianello quadrature nodes and weights for all polygons
* MESH.SHP --> shape functions on the ref. element
* MESH.SHG --> 4-D tensor of the shape functions gradient
%}


% Load Connectivity Matrix
for p = 1:mesh.NP
    mesh.conn(p,1:4) = mesh.polygon(p).vertices(:);  
end

% Load mesh dofs (vertices)
for v = 1:mesh.NV
mesh.xv(v) = mesh.vertex(v).x;
mesh.yv(v) = mesh.vertex(v).y;
end

% We use Vianello's quadrature rule  which integrate exactly polynomials of
% degree n
n = 2;


for ip = 1:mesh.NP
    
    % Storing the surface of each polygon
     mesh.area(ip,1) = mesh.polygon(ip).area;
    
    [xq,yq,wq]=quadrature_vianello(n,ip,mesh);
    
    % storing quadrature nodes and weights in the mesh structure
    mesh.polygon(ip).xq = xq;
    mesh.polygon(ip).yq = yq;
    mesh.polygon(ip).wq = wq;
   
    
    % Load here the shape functions and their gradient computed in the
    % quadrature nodes.
    % Notice: SHG(component,value of N in (xq,yq),N_i,Polygon)
      mesh.shp (:, 1, ip) = (1 - xq) .* (1 - yq);
      mesh.shp (:, 2, ip) = (xq) .* (1 - yq);
      mesh.shp (:, 3, ip) = (xq) .* (yq);
      mesh.shp (:, 4, ip) = (1 - xq) .* (yq);

      mesh.shg (1, :, 1, ip) = (-1+yq);
      mesh.shg (2, :, 1, ip) = (-1+xq);
      mesh.shg (1, :, 2, ip) = (1 - yq);
      mesh.shg (2, :, 2, ip) = -(xq);
      mesh.shg (1, :, 3, ip) = (yq);
      mesh.shg (2, :, 3, ip) = (xq);
      mesh.shg (1, :, 4, ip) = -(yq);
      mesh.shg (2, :, 4, ip) = (1 - xq);
    
end

return