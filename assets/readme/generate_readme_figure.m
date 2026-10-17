function generate_readme_figure
% Reproduce one cortical workflow figure using the bundled MATLAB functions.
% First run: python3 assets/readme/prepare_demo_surfaces.py
% Then in MATLAB, from the repository root:
%   addpath(fullfile(pwd,'assets','readme')); generate_readme_figure
% ASCII copies avoid the compressed GIfTI reader's Java dependency.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);
outdir = fullfile(root,'assets','readme');
builddir = fullfile(root,'build','readme');
p = gifti(fullfile(builddir,'pial.surfL.gii'));
q = gifti(fullfile(builddir,'fs_LR.32k.L.flat.surf.gii'));
s = struct('vertices',double(p.vertices),'faces',double(p.faces));
flat = struct('vertices',double(q.vertices),'faces',double(q.faces));
beta = surf2quadratic(s);
curv = surf2curvature(beta);
h = curv.mean(:);
[A,C] = FEM(s);
one = ones(size(h)); area = one'*A*one;
hc = h - (one'*A*h)/area;
lambda = [1 0.1 0.01 0];
U = zeros(numel(h),numel(lambda));
J = cell(1,numel(lambda));
checks.vertices = size(s.vertices,1);
checks.faces = size(s.faces,1);
checks.flat_vertices = size(flat.vertices,1);
checks.flat_faces = size(flat.faces,1);
checks.lambda = lambda;
checks.centered_source_mass_mean = (one'*A*hc)/area;
checks.stiffness_constant_residual = norm(C*one)/norm(C,'fro');
for k = 1:numel(lambda)
    [U(:,k),J{k}] = FEM_poisson_regualarized(s,h,lambda(k));
    if lambda(k)>0
        raw = (C+lambda(k)*A)\(A*hc);
    else
        b = A*one;
        sol = [C b;b' 0]\[A*hc;0];
        raw = sol(1:end-1);
    end
    gain = sqrt((hc'*A*hc)/area)/(sqrt((raw'*A*raw)/area)+eps);
    expected = gain*raw;
    expected = expected-(one'*A*expected)/area;
    checks.solver_agreement(k) = norm(U(:,k)-expected)/norm(expected);
    checks.scaled_system_residual(k) = norm((C+lambda(k)*A)*U(:,k)-gain*A*hc)/norm(gain*A*hc);
    checks.potential_mass_mean(k) = (one'*A*U(:,k))/area;
    checks.rms_ratio(k) = sqrt((U(:,k)'*A*U(:,k))/(hc'*A*hc));
end

% Independently check the returned vector orientation using edge differences.
% At this source version Jf points uphill; flux_to_display reverses Jf.
F=s.faces; V=s.vertices; u=U(:,2); Jf=J{2};
e1=V(F(:,2),:)-V(F(:,1),:); e2=V(F(:,3),:)-V(F(:,1),:);
d1=u(F(:,2))-u(F(:,1)); d2=u(F(:,3))-u(F(:,1));
checks.returned_vector_gradient_residual = norm([sum(Jf.*e1,2)-d1;sum(Jf.*e2,2)-d2])/norm([d1;d2]);
n=cross(e1,e2,2); n=n./sqrt(sum(n.^2,2));
checks.vector_tangency_residual = norm(sum(Jf.*n,2))/norm(Jf,'fro');
checks.display_limit = [-0.5 0.5];
checks.fraction_clipped = mean(abs([hc U(:,1:3)])>0.5,1);
assert(checks.vertices==checks.flat_vertices);
assert(max(checks.solver_agreement)<1e-9);
assert(max(checks.scaled_system_residual)<1e-8);
assert(max(abs(checks.rms_ratio-1))<1e-9);
assert(checks.returned_vector_gradient_residual<1e-9);
assert(checks.vector_tangency_residual<1e-9);

navy=[.10 .18 .29]; muted=[.35 .41 .49];
f=figure('Visible','off','Color','white','Position',[80 80 1800 1140]);
set(f,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial', ...
    'DefaultAxesColor','white','DefaultAxesXColor',navy,'DefaultAxesYColor',navy, ...
    'PaperPositionMode','auto','InvertHardcopy','off');
label(f,[.03 .940 .94 .045],'Poisson flow: from cortical geometry to a folding potential',27,navy);
label(f,[.03 .895 .94 .035],sprintf('Bundled left hemisphere  |  %s vertices  |  %s triangular faces  |  Computed in MATLAB', ...
    '32,492','64,980'),15,muted);
cols=[.03 .355 .68];
names={'A  Cortical surface','B  Mean-centered curvature','C  Potential and downhill flow', ...
       'D  More local folding detail','E  Broader folding organization','F  The same potential on a flat map'};
sub={'Pial mesh: the computational domain','Local quadratic surface fitting','Screening parameter \lambda = 0.1', ...
     '\lambda = 1','\lambda = 0.01','\lambda = 0.1; scalar values transferred by vertex'};
for k=1:6
    row=ceil(k/3); col=mod(k-1,3)+1;
    yy=.837-(row-1)*.373;
    label(f,[cols(col) yy .30 .04],names{k},19,navy);
    lab=label(f,[cols(col) yy-.035 .30 .032],sub{k},14,muted);
    set(lab,'Interpreter','tex');
end
maps={zeros(size(h)),hc,U(:,2),U(:,1),U(:,3),U(:,2)};
for k=1:6
    row=ceil(k/3); col=mod(k-1,3)+1;
    ax(k)=axes(f,'Position',[cols(col) .523-(row-1)*.373 .285 .266]);
    if k==6
        drawsurface(ax(k),flat,maps{k},true);
    else
        drawsurface(ax(k),s,maps{k},false);
    end
    if k==1
        colormap(ax(k),[.72 .75 .78]);
    end
    if k==3
        axes(ax(k));
        step=55; arrowScale=2.4;
        flux_to_display(s,Jf,step,arrowScale,'face');
        qu=findobj(ax(k),'Type','quiver');
        set(qu,'Color',[.04 .27 .17],'LineWidth',.7,'MaxHeadSize',1.3);
        displayed=[qu.UData(:) qu.VData(:) qu.WData(:)];
        ix=1:step:size(Jf,1);
        target=-arrowScale*Jf(ix,:)./(sqrt(sum(Jf(ix,:).^2,2))+eps);
        checks.display_arrow_agreement=norm(displayed-target,'fro')/norm(target,'fro');
        assert(checks.display_arrow_agreement<1e-12);
    end
end
cb=colorbar(ax(2),'southoutside');
cb.Position=[.355 .094 .285 .015]; cb.Color=navy; cb.FontSize=12;
cb.Ticks=[-.5 -.25 0 .25 .5];
cb.Label.String='Curvature / RMS-matched potential: shared display scale';
cb.Label.Color=navy; cb.Label.FontSize=12;
annotation(f,'line',[.03 .97],[.055 .055],'Color',[.83 .87 .90]);
label(f,[.03 .009 .94 .040], ...
    'All fields are computed on the pial surface. Arrows show direction only. The flat panel displays scalar values, not transported vectors. Example data; no group comparison.',12,muted);
print(f,fullfile(outdir,'poisson-workflow.png'),'-dpng','-r170');
close(f);
checks.matlab=version;
checks.source_commit='d144e3a720adb4c565b2c238c287e075f92ad0c3';
fid=fopen(fullfile(outdir,'numerical-checks.json'),'w');
fwrite(fid,jsonencode(checks,PrettyPrint=true)); fclose(fid);
save(fullfile(builddir,'figure-data.mat'),'s','flat','h','hc','lambda','U','J','checks');
disp(checks);
fprintf('Single cortical workflow figure exported successfully.\n');
end

function hLabel=label(f,pos,str,sz,color)
hLabel=annotation(f,'textbox',pos,'String',str,'EdgeColor','none','FontName','Arial', ...
    'FontSize',sz,'Color',color,'HorizontalAlignment','left', ...
    'VerticalAlignment','middle','Interpreter','none');
end

function drawsurface(ax,s,value,isflat)
axes(ax);
p=figure_trimesh(s,value,'rwb');
colormap(ax,diverging()); clim(ax,[-.5 .5]);
set(p,'AmbientStrength',.65,'DiffuseStrength',.55,'SpecularStrength',.06);
if isflat
    view(ax,2); lighting(ax,'none');
else
    view(ax,[-90 0]); camlight(ax,'headlight'); lighting(ax,'gouraud');
end
axis(ax,'equal'); axis(ax,'tight'); axis(ax,'off');
end

function c=diverging
t=linspace(0,1,128)';
left=[.14 .32 .68]; mid=[.98 .98 .97]; right=[.76 .18 .12];
c=[left.*(1-t)+mid.*t;mid.*(1-t)+right.*t];
end
