function h_mc = curvature_mean_centering(surf, h)
%function h_mc = curvature_mean_centering(surf, h)
%
%The function mean-centers a vertexwise curvature field h(v) on a closed
%triangulated surface with respect to local surface area a_v, i.e.,
%
%   h(v) <- h(v) - (sum_v a_v h(v)) / (sum_v a_v).
%
%This removes the constant (null-space) component and ensures solvability of
%downstream Poisson/diffusion equations on a closed manifold.
%
% INPUT
%  surf.vertices : n x 3 vertex coordinates
%  surf.faces    : m x 3 triangulation (1-based indices)
%  h             : n x 1 (or 1 x n) curvature field (e.g., curvatures.mean')
%
% OUTPUT
%  h_mc          : n x 1 mean-centered curvature field
%
% (C) Moo K. Chung mkchung@wisc.edu
% University of Wisconsin-Madison

V = surf.vertices;
F = surf.faces;

nvertex = size(V,1);

h = h(:);                                         % <— ensure column vector
% error checking for size mismatch should go here as comments.

% ---- barycentric (lumped) vertex areas a_v
cp   = cross(V(F(:,2),:) - V(F(:,1),:), ...
             V(F(:,3),:) - V(F(:,1),:), 2);       % <— added: explicit cross vector
Atri = 0.5 * sqrt(sum(cp.^2, 2));                 % <— fixed: unambiguous triangle area

a = zeros(nvertex,1);
a(F(:,1)) = a(F(:,1)) + Atri/3;
a(F(:,2)) = a(F(:,2)) + Atri/3;
a(F(:,3)) = a(F(:,3)) + Atri/3;

% ---- area-weighted mean and mean-centering
hbar = (a' * h) / (sum(a) + eps);                 % <— area-weighted mean
h_mc = h - hbar;                                  % <— mean-centered field
end