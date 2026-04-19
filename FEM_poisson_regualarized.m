function [u,Jf] = FEM_poisson_regualarized(surf,h,lambda)
%function [u,Jf] = FEM_poisson_regualarized(surf,h,lambda)
%
%Solves (screened) Poisson on a closed triangle mesh and returns face-wise flux:
%   (C + lambda*A) * u = A * h
%   J = -grad(u)  (computed per face)
%
%NUMERICAL NOTES (closed surfaces):
%   1) For lambda > 0, (C + lambda*A) is SPD and u is unique (no constraint needed).
%   2) For lambda = 0, C is singular; we fix gauge by enforcing 1'*A*u = 0.
%   3) To keep amplitude comparable across subjects/meshes and across lambda,
%      we normalize u to have the same area-weighted RMS as h.                 % <— added
%
%INPUT
% surf   : surf.vertices (n x 3), surf.faces (m x 3)
% h      : n x 1 vertexwise source/sink field (e.g., mean curvature)
% lambda : scalar >= 0 (lambda=0 gives standard Poisson; lambda>0 localizes)
%
%OUTPUT
% u  : n x 1 potential (gauge-fixed if lambda=0; unique if lambda>0)
% Jf : m x 3 face-wise flux vectors (points toward sinks)
%
% (C) Moo K. Chung mkchung@wisc.edu
% University of Wisconsin-Madison

if nargin < 3
    lambda = 0.1;
end

[A,C] = FEM(surf);

h   = h(:);
one = ones(size(h));

% ---- enforce zero-mean constraint
hbar = (one' * A * h) / (one' * A * one);
h    = h - hbar;                                                     % <— keep

b = A*h;

if lambda > 0
    % ---- regularized solve (unique; no constraint needed)
    K0 = C + lambda*A;
    u  = K0 \ b;

else
    % ---- standard Poisson: zero-mean constraint needed
    g   = A*one; 
    K   = [C, g;
           g', 0];
    rhs = [b; 0];
    sol = K \ rhs;
    u   = sol(1:end-1);
end


% Gauge condition: scale normalization to keep u comparable to original data scale
% Match area-weighted RMS(u) to area-weighted RMS(h).                           % <— added
h_rms = sqrt( (h' * A * h) / (one' * A * one) );
u_rms = sqrt( (u' * A * u) / (one' * A * one) );
u     = u * (h_rms / (u_rms + eps));                                           % <— added

% ---- Remove numerical mean drift of u (does not change grad(u))
ubar = (one' * A * u) / (one' * A * one);
u    = u - ubar;                                                               % <— added

% ---- face-wise flux J = -grad(u)
V = surf.vertices;
F = surf.faces;

i = F(:,1); j = F(:,2); k = F(:,3);
Pi = V(i,:); Pj = V(j,:); Pk = V(k,:);
ui = u(i);   uj = u(j);   uk = u(k);

e1 = Pj - Pi;
e2 = Pk - Pi;
n  = cross(e1,e2,2);
A2 = sqrt(sum(n.^2,2));                    % 2*area
n_unit = n ./ A2;

grad_phi_i = cross(n_unit, (Pj-Pk), 2) ./ A2;
grad_phi_j = cross(n_unit, (Pk-Pi), 2) ./ A2;
grad_phi_k = cross(n_unit, (Pi-Pj), 2) ./ A2;

grad_u = grad_phi_i .* ui + grad_phi_j .* uj + grad_phi_k .* uk;
Jf     = -grad_u;

end