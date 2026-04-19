function flux_to_display(surf,Jf,step,scale,mode)
%function flux_to_display(surf,Jf,step,scale,mode)
%
%Dense flux visualization OVERLAY on existing surface plot.
%Supports face-centroid arrows (default) or vertex-based arrows (denser).
%
% INPUT
%  surf.vertices : n x 3
%  surf.faces    : m x 3
%  Jf            : m x 3 face-wise flux vectors
%  step          : subsampling step (default 5; smaller => denser)
%  scale         : arrow length scale (recommended 1~2; default 1.5)
%  mode          : 'face' or 'vertex' (default 'vertex')
%
% OUTPUT
%  draws arrows on current axes
%
% (C) Moo K. Chung mkchung@wisc.edu
% University of Wisconsin-Madison

if nargin < 3, step  = 5;    end
if nargin < 4, scale = 1.5;  end
if nargin < 5, mode  = 'vertex'; end

V = surf.vertices;
F = surf.faces;

hold on;

switch lower(mode)
case 'face'

    % ---- face centroids
    Cc = (V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:))/3;

    idx = 1:step:size(F,1);
    P   = Cc(idx,:);
    J   = Jf(idx,:);

case 'vertex'

    % ---- convert face-wise flux to vertex-wise flux by accumulation
    n = size(V,1);
    Jv  = zeros(n,3);
    cnt = zeros(n,1);

    for k = 1:3
        id = F(:,k);
        Jv(id,:) = Jv(id,:) + Jf;
        cnt(id)  = cnt(id) + 1;
    end
    Jv = Jv ./ max(cnt,1);

    idx = 1:step:size(V,1);
    P   = V(idx,:);
    J   = Jv(idx,:);

otherwise
    % mode must be 'face' or 'vertex'  % <— comment only
    P = []; J = [];
end

% ---- direction-only arrows (unit vectors) so scale is meaningful (1~2)
% ---- enforce flow from high to low (flux = -grad u)
J  = -J;                                   
Jn = sqrt(sum(J.^2,2)) + eps;
J  = J ./ Jn;

% ---- draw (do NOT reset colormap / do NOT wash out surface)
quiver3(P(:,1),P(:,2),P(:,3), ...
        scale*J(:,1),scale*J(:,2),scale*J(:,3), ...
        0,'g','LineWidth',1,'MaxHeadSize',2); 

end