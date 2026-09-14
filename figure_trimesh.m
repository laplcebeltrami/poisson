function h = figure_trimesh(surf,value,cmap)
%
% figure_trimesh(surf,value,cmap)
% Color mesh vertices of "surf" using given "value" with colormap "cmap".
%
% cmap options:  'rwb','ryw','mwg','rywb','rywgb','rw','oyw','w','default'
%
% (C) 2008. Moo K. Chung, 
%     University of Wisconsin-Madison
% 
%  Contact: email: mkchung@wisc.edu
%
% Update history: Apr. 29, 2008  created
%                 Feb. 02, 2009  additional color scheme by S.G. Kim 
%                 Nov. 26, 2018  document updated
%                 May 1, 2019  document updated
%                 Jan. 7, 2026 simplified colormap construction (no hard-coded RGB tables)

if nargin < 3, cmap = 'default'; end

coord = surf.vertices;
tri   = surf.faces;

hold on;
% ---- plot
h=trisurf(tri,coord(:,1),coord(:,2),coord(:,3),value,'EdgeColor','none');

set(h,'UserData',value); %<--new

daspect([1 1 1]);  axis tight off;
lighting gouraud; material shiny; shading interp;
set(gcf,'Color','w'); colorbar off;

% ---- colormap
cm = get_cmap(cmap);                                % <— added: force numeric Nx3 RGB map
colormap(gca, cm);                                  % <— FIX: apply to CURRENT AXES, not globally


% =========================
% Local colormap factory
% =========================
function cm = get_cmap(name)
name = lower(string(name));

switch name
    case "default"
        cm = colormap;  % whatever MATLAB currently has

    case "w"
        cm = [1 1 1];

    case "rw"
        cm = ramp([1 1 1],[1 0 0],100);                 % white -> red

    case "rwb"
        cm = diverge([0 0 1],[1 1 1],[1 0 0],101);       % blue -> white -> red

    case "ryw"
        cm = diverge([1 0 0],[1 1 0],[1 1 1],101);       % red -> yellow -> white

    case 'rywb'
    % deep blue -> blue -> white -> yellow -> red  (matches your rwb table)
    cm = ramp5([0 0 0.5625],[0 0 1],[1 1 1],[1 1 0],[1 0 0],129);  % <— FIX
    
    %case "rywb"
    %    cm = ramp4([0 0 0.5625],[0 0 1],[1 1 1],[1 0 0],129); % deep blue -> blue -> white -> red
    
   % case 'rywgb'
   %  % red -> yellow -> white (wide) -> green -> blue
   %  cm = [ ...
   %      ramp([1 0 0],[1 1 0],30); ...      % red -> yellow
   %      ramp([1 1 0],[1 1 1],25); ...      % yellow -> white
   %      repmat([1 1 1],20,1); ...          % <— WIDE WHITE BAND
   %      ramp([1 1 1],[0 1 0],25); ...      % white -> green
   %      ramp([0 1 0],[0 0 1],29) ];        % green -> blue

    case 'rywgb'
    % blue -> green -> white -> yellow -> red
    cm = [ ...
        ramp([0 0 1],[0 1 0],30); ...      % blue -> green
        ramp([0 1 0],[1 1 1],25); ...      % green -> white
        repmat([1 1 1],20,1); ...          % wide white band
        ramp([1 1 1],[1 1 0],25); ...      % white -> yellow
        ramp([1 1 0],[1 0 0],29) ];        % yellow -> red

    case 'oldrywb'

        cm =[0         0    0.5625
         0         0    0.6500
         0         0    0.7375
         0         0    0.8250
         0         0    0.9125
         0         0    1.0000
         0    0.0909    1.0000
         0    0.1818    1.0000
         0    0.2727    1.0000
         0    0.3636    1.0000
         0    0.4545    1.0000
         0    0.5455    1.0000
         0    0.6364    1.0000
         0    0.7273    1.0000
         0    0.8182    1.0000
         0    0.9091    1.0000
         0    1.0000    1.0000
    0.0859    0.9971    0.9954
    0.1718    0.9943    0.9907
    0.2578    0.9914    0.9861
    0.3437    0.9886    0.9815
    0.4296    0.9857    0.9768
    0.5155    0.9829    0.9722
    0.6014    0.9800    0.9676
    0.6873    0.9772    0.9629
    0.7733    0.9743    0.9583
    0.8592    0.9715    0.9537
    0.9451    0.9686    0.9490
    0.9451    0.9686    0.9490
    0.9451    0.9686    0.9490
    0.9451    0.9686    0.9490
    0.9451    0.9686    0.9490
    0.9451    0.9686    0.9490
    0.9451    0.9686    0.9490
    0.9451    0.9686    0.9490
    0.9451    0.9686    0.9490
    0.9497    0.9712    0.8699
    0.9542    0.9739    0.7908
    0.9588    0.9765    0.7118
    0.9634    0.9791    0.6327
    0.9680    0.9817    0.5536
    0.9725    0.9843    0.4745
    0.9771    0.9869    0.3954
    0.9817    0.9895    0.3163
    0.9863    0.9922    0.2373
    0.9908    0.9948    0.1582
    0.9954    0.9974    0.0791
    1.0000    1.0000         0
    1.0000    0.9091         0
    1.0000    0.8182         0
    1.0000    0.7273         0
    1.0000    0.6364         0
    1.0000    0.5455         0
    1.0000    0.4545         0
    1.0000    0.3636         0
    1.0000    0.2727         0
    1.0000    0.1818         0
    1.0000    0.0909         0
    1.0000         0         0
    0.9000         0         0
    0.8000         0         0
    0.7000         0         0
    0.6000         0         0
    0.5000         0         0];


    case "mwg"
        cm = diverge([0.8 0 0.8],[1 1 1],[0 0.6 0.6],101);     % magenta -> white -> green

    case "oyw"
        % your "orange-yellow-white" idea: keep ryw structure but push last part to orange
        base = diverge([1 0 0],[1 1 0],[1 1 1],101);           % red->yellow->white
        startIndex = 37;
        orangeTarget = [1 0.5 0];
        cm = base;
        cm(startIndex:end,:) = ramp(cm(startIndex,:), orangeTarget, size(cm,1)-startIndex+1);

    otherwise
        cm = colormap;
end


% =========================
% Helpers
% =========================

function cm = diverge(cL,cM,cR,n)
% n should be odd so that exact center exists
if mod(n,2)==0, n = n+1; end
n2 = (n+1)/2;
cm = [ramp(cL,cM,n2); ramp(cM,cR,n2)];
cm(n2,:) = []; % remove duplicated middle row

function cm = ramp(c1,c2,n)
t  = linspace(0,1,n)';
cm = (1-t).*c1 + t.*c2;


function cm = ramp4(c1,c2,c3,c4,n)
% piecewise linear through 4 anchors
n1 = floor(n/3); 
n2 = floor(n/3); 
n3 = n - n1 - n2;
cm = [ramp(c1,c2,n1); ramp(c2,c3,n2); ramp(c3,c4,n3)];


function cm = ramp5(c1,c2,c3,c4,c5,n)
% piecewise linear interpolation through 5 anchor colors
n1 = floor(n/4);
n2 = floor(n/4);
n3 = floor(n/4);
n4 = n - n1 - n2 - n3;

cm = [ ramp(c1,c2,n1); ...
       ramp(c2,c3,n2); ...
       ramp(c3,c4,n3); ...
       ramp(c4,c5,n4) ];

% function figure_trimesh(surf,value, cmap)
% %
% % function figure_trimesh(surf,value, cmap)
% % 
% % Color mesh vertices of "surf" using given "value" using colormap "cmap". 
% % If "cmap" is not given, the MATLAB default colormap will be used.
% %
% % cmap options:  'rwb' : red-white-blue
% %                'ryw' : red-yellow-white
% %                'mwg' : magenta-white-green added by SG Kim in 2009.
% %                'rywb': red-yellow-white-blue
% %
% %
% % (C) 2008. Moo K. Chung, S.G. Kim 
% %     University of Wisconsin-Madison
% % 
% %  Contact: email: mkchung@wisc.edu
% %
% % 
% % Update history: Apr. 29, 2008  created
% %                 Feb. 02, 2009  additional color scheme added
% %                 Nov. 26, 2018  document updated
% %                 May 1, 2019  document updated
% 
% 
% %red-white-blue
% rwb=[ 0         0    1.0000
%     0.0512    0.0565    1.0000
%     0.1024    0.1131    1.0000
%     0.1536    0.1696    1.0000
%     0.2048    0.2261    1.0000
%     0.2560    0.2826    1.0000
%     0.3072    0.3392    1.0000
%     0.3584    0.3957    1.0000
%     0.4096    0.4522    1.0000
%     0.4608    0.5087    1.0000
%     0.5120    0.5653    1.0000
%     0.5632    0.6218    1.0000
%     0.6145    0.6783    1.0000
%     0.6657    0.7348    1.0000
%     0.7169    0.7914    1.0000
%     0.7326    0.8030    1.0000
%     0.7483    0.8146    1.0000
%     0.7641    0.8261    1.0000
%     0.7798    0.8377    1.0000
%     0.7955    0.8493    1.0000
%     0.8112    0.8609    1.0000
%     0.8270    0.8725    1.0000
%     0.8427    0.8841    1.0000
%     0.8584    0.8957    1.0000
%     0.8742    0.9073    1.0000
%     0.8899    0.9189    1.0000
%     0.9056    0.9305    1.0000
%     0.9214    0.9420    1.0000
%     0.9371    0.9536    1.0000
%     0.9528    0.9652    1.0000
%     0.9685    0.9768    1.0000
%     0.9843    0.9884    1.0000
%     1.0000    1.0000    1.0000
%     1.0000    1.0000    0.9333
%     1.0000    1.0000    0.8667
%     1.0000    1.0000    0.8000
%     1.0000    1.0000    0.7333
%     1.0000    1.0000    0.6667
%     1.0000    1.0000    0.6000
%     1.0000    1.0000    0.5333
%     1.0000    1.0000    0.4667
%     1.0000    1.0000    0.4000
%     1.0000    1.0000    0.3333
%     1.0000    1.0000    0.2667
%     1.0000    1.0000    0.2000
%     1.0000    1.0000    0.1333
%     1.0000    1.0000    0.0667
%     1.0000    1.0000         0
%     0.9904    0.9475         0
%     0.9809    0.8951         0
%     0.9713    0.8426         0
%     0.9618    0.7902         0
%     0.9522    0.7377         0
%     0.9426    0.6853         0
%     0.9331    0.6328         0
%     0.9235    0.5804         0
%     0.9140    0.5279         0
%     0.9044    0.4755         0
%     0.8949    0.4230         0
%     0.8853    0.3706         0
%     0.8757    0.3181         0
%     0.8662    0.2657         0
%     0.8566    0.2132         0
%     0.8471    0.1608         0];
% 
% %red-yellow-white
% ryw =[ 1.0000    1.0000    1.0000
%     1.0000    0.9991    0.9977
%     1.0000    0.9982    0.9954
%     1.0000    0.9972    0.9931
%     1.0000    0.9963    0.9908
%     1.0000    0.9954    0.9885
%     1.0000    0.9945    0.9862
%     1.0000    0.9935    0.9839
%     1.0000    0.9926    0.9815
%     1.0000    0.9917    0.9792
%     1.0000    0.9908    0.9769
%     1.0000    0.9899    0.9746
%     1.0000    0.9889    0.9723
%     1.0000    0.9880    0.9700
%     1.0000    0.9871    0.9677
%     1.0000    0.9862    0.9654
%     1.0000    0.9852    0.9631
%     1.0000    0.9843    0.9608
%     1.0000    0.9834    0.9585
%     1.0000    0.9825    0.9562
%     1.0000    0.9815    0.9539
%     1.0000    0.9806    0.9516
%     1.0000    0.9797    0.9493
%     1.0000    0.9788    0.9469
%     1.0000    0.9779    0.9446
%     1.0000    0.9769    0.9423
%     1.0000    0.9760    0.9400
%     1.0000    0.9751    0.9377
%     1.0000    0.9742    0.9354
%     1.0000    0.9732    0.9331
%     1.0000    0.9723    0.9308
%     1.0000    0.9714    0.9285
%     1.0000    0.9705    0.9262
%     1.0000    0.9696    0.9239
%     1.0000    0.9686    0.9216
%     1.0000    0.9709    0.8557
%     1.0000    0.9731    0.7899
%     1.0000    0.9754    0.7241
%     1.0000    0.9776    0.6583
%     1.0000    0.9798    0.5924
%     1.0000    0.9821    0.5266
%     1.0000    0.9843    0.4608
%     1.0000    0.9866    0.3950
%     1.0000    0.9888    0.3291
%     1.0000    0.9910    0.2633
%     1.0000    0.9933    0.1975
%     1.0000    0.9955    0.1317
%     1.0000    0.9978    0.0658
%     1.0000    1.0000         0
%     0.9861    0.9237         0
%     0.9722    0.8474         0
%     0.9583    0.7711         0
%     0.9444    0.6948         0
%     0.9305    0.6185         0
%     0.9166    0.5422         0
%     0.9027    0.4660         0
%     0.8888    0.3897         0
%     0.8749    0.3134         0
%     0.8610    0.2371         0
%     0.8471    0.1608         0
%     0.7686    0.1480    0.0020
%     0.6902    0.1353    0.0039
%     0.6118    0.1225    0.0059
%     0.5333    0.1098    0.0078];
% 
% 


% 
% % magenta-white-green
% mwg =[   0    0.6000    0.6000
%     0.0034    0.6347    0.6347
%     0.0067    0.6695    0.6695
%     0.0101    0.7042    0.7042
%     0.0135    0.7389    0.7389
%     0.0168    0.7737    0.7737
%     0.0202    0.8084    0.8084
%     0.0235    0.8431    0.8431
%     0.0716    0.8520    0.8520
%     0.1197    0.8609    0.8609
%     0.1678    0.8698    0.8698
%     0.2159    0.8787    0.8787
%     0.2641    0.8876    0.8876
%     0.3122    0.8965    0.8965
%     0.3603    0.9054    0.9054
%     0.4084    0.9142    0.9142
%     0.4565    0.9231    0.9231
%     0.5046    0.9320    0.9320
%     0.5527    0.9409    0.9409
%     0.6008    0.9498    0.9498
%     0.6489    0.9587    0.9587
%     0.6970    0.9676    0.9676
%     0.7451    0.9765    0.9765
%     0.7706    0.9788    0.9788
%     0.7961    0.9812    0.9812
%     0.8216    0.9835    0.9835
%     0.8471    0.9859    0.9859
%     0.8725    0.9882    0.9882
%     0.8980    0.9906    0.9906
%     0.9235    0.9929    0.9929
%     0.9490    0.9953    0.9953
%     0.9745    0.9976    0.9976
%     1.0000    1.0000    1.0000
%     0.9980    0.9784    0.9980
%     0.9961    0.9569    0.9961
%     0.9941    0.9353    0.9941
%     0.9922    0.9137    0.9922
%     0.9902    0.8922    0.9902
%     0.9882    0.8706    0.9882
%     0.9863    0.8490    0.9863
%     0.9843    0.8275    0.9843
%     0.9824    0.8059    0.9824
%     0.9804    0.7843    0.9804
%     0.9731    0.7349    0.9731
%     0.9658    0.6855    0.9658
%     0.9584    0.6361    0.9584
%     0.9511    0.5867    0.9511
%     0.9438    0.5373    0.9438
%     0.9365    0.4878    0.9365
%     0.9292    0.4384    0.9292
%     0.9218    0.3890    0.9218
%     0.9145    0.3396    0.9145
%     0.9072    0.2902    0.9072
%     0.8999    0.2408    0.8999
%     0.8925    0.1914    0.8925
%     0.8852    0.1420    0.8852
%     0.8779    0.0925    0.8779
%     0.8706    0.0431    0.8706
%     0.8588    0.0360    0.8588
%     0.8471    0.0288    0.8471
%     0.8353    0.0216    0.8353
%     0.8235    0.0144    0.8235
%     0.8118    0.0072    0.8118
%     0.8000         0    0.8000];
% 
% 
% 
% % Modify the later part of the colormap to transition to orange
% % Orange is typically represented as RGB [1, 0.5, 0], so we target that.
% orangeTarget = [1, 0.5, 0];
% 
% % Index where the modification starts (this can be adjusted)
% startIndex = 37; % Index from where orange color starts
% 
% % Number of steps from startIndex to the end of the colormap
% nSteps = size(ryw, 1) - startIndex + 1;
% 
% oyw = ryw;
% % Create a linear transition from the color at startIndex to orange
% for i = startIndex:size(ryw, 1)
%     fraction = (i - startIndex) / (nSteps - 1);
%     oyw(i, :) = (1 - fraction) * oyw(startIndex, :) + fraction * orangeTarget;
% end
% 
% rw = red_to_white(100);
% 
% 
% if nargin < 3
%     cmap='default';
% end
% 
% edge_c = 'k';
% 
% switch cmap
%     case 'default'
%         color_map=cmap;
%     case 'rw'
%         color_map=rw;
%     case 'rwb'
%         color_map=rwb;
%     case 'ryw'
%         color_map=ryw;
%     case 'rywb'
%         color_map=rywb;
%     case 'mwg'
%         color_map=mwg;
%     case 'oyw'
%         color_map=oyw;
%     case 'w'
%         color_map = [1 1 1]; % pure white color  
%         edge_c = 'w';
% end;
% 
% 
% 
% %edge_c = 'k';
% edge_c = 'none';
% 
% coord=surf.vertices;
% tri=surf.faces;
% 
% trisurf(tri,coord(:,1), coord(:,2),coord(:,3),value, 'EdgeColor', edge_c);colorbar;
% 
% 
% colormap(color_map);
% 
% daspect([1 1 1]);
% %view(3); 
% axis tight;
% axis vis3d off;
% 
% 
% 
% 
% %camlight headlight; 
% lighting gouraud %phong %gouraud; 
% material shiny %metal %shiny 
% shading interp;
% %camlight (-120,15);
% %alpha(1)
% 
% set(gcf,'Color','w')
% 
% colorbar('off')
% 
% %------------
% function cmap = red_to_white(steps)
%     cmap = zeros(steps, 3);  % Initialize the colormap
%     for i = 1:steps
%         t = (i - 1) / (steps - 1);  % Normalized position in the gradient
%         cmap(steps-i+1, :) = [1, t, t];  % Linear interpolation between red and white
%     end
% 
% 
% 
