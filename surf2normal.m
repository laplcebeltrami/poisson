function surfnormal =surf2normal(surf);
%
% normal =surf2normal(surf);
%
% The function computes the unit normal vector of triangulated surface at vertex.
%
%
% The algorithm is based on 
%
% [1] Chung, M.K., 2013 Computational Neuroanatomy: The Methods, World
% Scientific Publishing, pages 100-109
%
% 
% [2] Chung. M.K. 2001, Statistical Morphometry in Computional Neuroanatomy, 
%     PhD Thesis, McGill University, Montreal. 
%     http://www.stat.wisc.edu/~mchung/papers/thesis.pdf
%
%
% (C) 2019 Moo K. Chung
% University of Wisconsin-Madiosn
%
% mkchung@wisc.edu
%
% 
% Update history
% 2019 August 23: Modified from Getarea.m
% 2019 August 26: Gives unit normal vector

area = faceArea(surf); %area of face
normals = faceNormal(surf); %normal of face

faces = surf.faces;
vertices = surf.vertices;
nvertices = size(vertices,1);

surfnormal=zeros(nvertices,3);
for i=1:nvertices
   [ind_tri, ind_ver] = find(faces ==i); 
   areas = area(ind_tri); %areas as weight
   
   %Old code gives normal vector but will not give unit normal vector.
   %areas = areas/sum(areas); %normalize areas
   %surfnormal(i,:) = sum(areas.*normals(ind_tri,:),1); 
   
   unitnormal = sum(areas.*normals(ind_tri,:),1);
   unitnormal = unitnormal/norm(unitnormal);
   surfnormal(i,:) = unitnormal; 
end





function area = faceArea(surf)
%function area = faceArea(surf)
%
%The function computes the area of each triangle face.
%
%(C) Moo K. Chung
% University of Wisconsin-Madison
%
%2019 August 23


faces=surf.faces;
coord=surf.vertices;

v1 = coord(faces(:,2),:) - coord(faces(:,1),:); %first edges
%v1 = v1./vecnorm(v1,2,2);  %unit vector 
v2 = coord(faces(:,3),:) - coord(faces(:,1),:); %second edges
%v2 = v2./vecnorm(v2,2,2);  %unit vector
%angels = acos(sum(v1.*v2,2)); %cosine angle between v1 and v2 

area = vecnorm(cross(v1,v2),2,2);



function normals = faceNormal(surf)
% function normals = faceNormal(surf)
%
% This function calculates the normal vectors for each face of a 2D surface.
% The function calculates the vectors of the first and second edges for each face.
% It then calculates the cross product of these two vectors to get the normal vector 
% for each face. Finally, it normalizes each normal vector.
%
% Input:
% surf - a structure with two fields: 'vertices' and 'faces'. 'vertices' is a 
%        Nx3 matrix that contains the coordinates (x, y, z) of each vertex, where
%        N is the number of vertices. 'faces' is a Mx3 matrix that contains the 
%        indices of the vertices that make up each face, where M is the number of faces.
%
% Output:
% normals - a Mx3 matrix that contains the normal vector (nx, ny, nz) for each face. 
%
%(C) Moo K. Chung
% University of Wisconsin-Madison
%2019 August 23

nodes=surf.vertices; 
faces=surf.faces; 

coord=surf.vertices; 

faces = surf.faces;

v1 = coord(faces(:,2),:) - coord(faces(:,1),:); %first edges 
v2 = coord(faces(:,3),:) - coord(faces(:,1),:); %second edges 

%Calculate the cross product of the two vectors (which gives the normal vector), and normalize it.
normals = cross(v1, v2); 
normals = normals./vecnorm(normals')';


