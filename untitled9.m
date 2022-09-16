close all; clear
global HeaterTemp Heater Bubblef HeaterText Oil v3 v4
figure('Color','w', 'Units', 'normalized', 'Position', [0.1300 0.1100 0.7750 0.8150], ....
    'keypressfcn',@fh_kpfcn);
axis([-4,10,-5,8]); daspect([1,1,1]); grid on; hold on;
%% Distillation Column
e = pi/40;
t1 = linspace(e, pi/2-e, 10); c1 = cos(t1); s1 = sin(t1);
t2 = linspace(e+pi/2, pi, 10); c2 = cos(t2); s2 = sin(t2);
t3 = linspace(pi, 1.5*pi-e, 10); c3 = cos(t3); s3 = sin(t3);
t4 = linspace(1.5*pi+e, 2*pi-e, 10); c4 = cos(t4); s4 = sin(t4);
v1 = 0*c1; v2 = 0*c2; v3 = 0*c3; v4 = 0*c4;
Tower.Fill = fill([c1,c2,c3,c4,c1(1)],[s1+5,s2+5,s3,s4,s1(1)+5],...
              [v1,v2,v3,v4,0],'EdgeAlpha',0, 'FaceAlpha',0.5);
Tower.Bound = plot([c1, nan, c2, c3, nan, c4, nan, c4(end),c1(1)],...
    [s1+5, nan, s2+5, s3, nan, s4, nan, sin(e),5-sin(e),],'k','LineWidth',5);
Tower.VelocityBias = @(bubble, a) VeloctyBiasY(bubble, a, 0, 5, 1.1);
Tower.Name = 'Tower';
%% Boiler
xb = 3+[c1,c2,c3,c4,c1(1)];
yb = [s1,s2,s3-1,s4-1,s1(1)]-1.4;
h = -1.2;
X = boilerlevel(xb, yb, h);
xs = linspace(X(1), X(2), 20);
ys = h+xs*0;
Boiler.Fill = plot([c1,c2,c3,c4,c1(1)]+3, [s1+1,s2+1,s3,s4,s1(1)+1]-2.4, ...
              'Color','none');
Boiler.Bound = plot([c1, nan, c2, nan, -1, c3, c4, c1(1)]+3, ...
    [s1, nan, s2, nan, -2*sin(e), s3-1, s4-1, s1(1)]-1.4, 'k','LineWidth',5);
Boiler.VelocityBias = @(bubble, a) VeloctyBias2(bubble, a);
Boiler.Name = 'Boiler';
Oil = fill([xs,3+[c3,c4]],[ys,[s3-1,s4-1]-1.4],...
            'b','EdgeAlpha',0, 'FaceAlpha',0.5);

%% Heater
t = pi/4*linspace(2,30,1001); 
OffTemp = 20; MaxTemp = 2000;
OffColor = [0,0,0]; MaxColor = [1,0,0];
xh = 3.2-(1-3*cos(t))/6; yh = 0.04*(t - 3*sin(t))-2.5;
xh = [5, xh, 5]; yh = [yh(1), yh, yh(end)];
Heater = plot(xh, yh,'k','LineWidth',2);

t = linspace(0, pi/2, 20); c = cos(t); s = sin(t);
xp1 = [3.1, 2.7 + 0.4*c, 0.8]; xp2 = [2.9, 2.7+0.2*c, 0.8];
yp1 = [-0.6,0.4*s-0.3, 0.1]; yp2 = [-0.6,0.2*s-0.3,-0.1];
Pipe1.Bound = plot([xp1, nan, flip(xp2), nan, xp1(1)],...
                  [yp1, nan, flip(yp2), nan, yp1(1)],...
                  'k', 'LineWidth',3);
Pipe1.Fill = fill([xp1, flip(xp2), xp1(1)], ...
    [yp1, flip(yp2), yp1(1)], 'w', 'EdgeAlpha',0, 'FaceAlpha',0.5);
Pipe1.VelocityBias = @(bubble, a) VeloctyBiasX(bubble, a, 0.8, 2.7, -1);
Pipe1.Name = 'Pipe1';

xp1 = [-0.1, 0.3 - 0.4*c, 5]; xp2 = [0.1, 0.3 - 0.2*c, 5];
yp1 = [5.8,0.4*s + 6.1, 6.5]; yp2 = [5.8,0.2*s + 6.1, 6.3];
Pipe2.Bound = plot([xp1, nan, flip(xp2), nan, xp1(1)],...
                  [yp1, nan, flip(yp2), nan, yp1(1)],...
                  'k', 'LineWidth',3);
Pipe2.Fill = fill([xp1, flip(xp2), xp1(1)], ...
    [yp1, flip(yp2), yp1(1)], 'w', 'EdgeAlpha',0, 'FaceAlpha',0.5);
Pipe2.VelocityBias = @(bubble, a) VeloctyBiasX(bubble, a, 0.3, 5, 1);
Pipe2.Name = 'Pipe2';

xp1 = [-0.1, 0.3 - 0.4*c, 2.2]; xp2 = [0.1, 0.3 - 0.2*c, 2.2];
yp1 = [-0.8,-1.2 - 0.4*s, -1.6]; yp2 = [-0.8,-1.2- 0.2*s, -1.4];
Pipe3.Bound = plot([xp1, nan, flip(xp2), nan, xp1(1)],...
                  [yp1, nan, flip(yp2), nan, yp1(1)],...
                  'k', 'LineWidth',3);
Pipe3.Fill = fill([xp1, flip(xp2), xp1(1)], ...
    [yp1, flip(yp2), yp1(1)], 'w', 'EdgeAlpha',0, 'FaceAlpha',0.5);
Pipe3.VelocityBias = @(bubble, a) VeloctyBiasX(bubble, a, 0.3, 2.2, 1);
Pipe3.Name = 'Pipe3';

xp1 = [0.9, 5]; xp2 = [0.9, 5];
yp1 = [4.9,4.9]; yp2 = [5.1,5.1];
Pipe4.Bound = plot([xp1, nan, flip(xp2), nan, xp1(1)],...
                  [yp1, nan, flip(yp2), nan, yp1(1)],...
                  'k', 'LineWidth',3);
Pipe4.Fill = fill([xp1, flip(xp2), xp1(1)], ...
    [yp1, flip(yp2), yp1(1)], 'w', 'EdgeAlpha',0, 'FaceAlpha',0.5);
Pipe4.VelocityBias = @(bubble, a) VeloctyBiasX(bubble, a, 1, 5, 1);
Pipe4.Name = 'Pipe4';


Components = [Pipe1, Pipe2, Pipe3, Pipe4, Tower, Boiler];
drawnow; Bubbles = []; Bubblef = 0; HeaterTemp = 20;
HeaterText = text(3,-3,['$',num2str(HeaterTemp),'^\circ C $'],...
            'interpreter', 'latex', 'Color','w', ...
            'HorizontalAlignment','center', 'FontSize',15);
t = 0; dt = 0.02;
vidfile = VideoWriter('Distiller6.avi');
open(vidfile);
for n = 1:2000
%while(1)
    t = t + dt;
    ys = Wavemaker(xs, h, t, Bubblef);
    Oil.YData = [ys,[s3-1,s4-1]-1.4];
    Tower.Fill.CData = [v1,v2,v3,v4,0];
    Bubbles = UpdateBubbles(Bubbles, xs, ys, dt, Components);
    if(rand < Bubblef)
        bx = 3 + 1.4*(rand-0.5); by = -2 + 0.8*(rand-0.5);
        bu = 0.5*(rand-0.5); bv = -0.5*rand; br = 0.05*rand;
        bubble = BubbleMaker(bx, by, bu, bv, br, xs, ys);
        Bubbles = [Bubbles, bubble];
    end
    drawnow
    frame = getframe(gcf);
    writeVideo(vidfile,frame);
end
close(vidfile);

function Bubbles = UpdateBubbles(Bubbles, xs, ys, dt, plthandles)
    I = zeros(numel(Bubbles),1);
    for i = 1:numel(Bubbles)
        bubble = Bubbles(i);
        if(bubble.position(1)<-2 ||bubble.position(1) >5 ...
           || bubble.position(2) < -5 || bubble.position(2)>8)
            delete(bubble.bob);
            I(i) = 1;
        elseif(bubble.collected)
            delete(bubble.bob);
            I(i) = 1;
        else
            bubble = bubble.Update(xs, ys, dt, plthandles);
            Bubbles(i) = bubble;
        end
    end
    Bubbles(I==1) = [];
end

function ys = Wavemaker(xs, h, t, Bubblef)
    ys = h + Bubblef*(0.2*sin(5*Bubblef*(7*xs - 3*t)) - 0.05*rand(size(xs)));
end

function [bl, p0] = Intersect(a, b, x, y)
    den0 = (b(1) - a(1)) * (x(2) - y(2)) - (b(2) - a(2)) * (x(1) - y(1));
    num1 = (x(1) - y(1)) * (a(2) - x(2)) - (x(2) - y(2)) * (a(1) - x(1));
    num2 = (a(1) - x(1)) * (b(2) - a(2)) - (a(2) - x(2)) * (b(1) - a(1));
    num1 = num1*(abs(num1) > 1e-10); num2 = num2*(abs(num2) > 1e-10);
    s = num1 / den0; t = num2 / den0; 
    bl = (0 < s && s < 1) && (0 < t && t < 1);
    p0 = a+s*(b-a);
end

function X = boilerlevel(xb, yb, h)
    p1 = [1.5,h]; p2 = [4.5,h]; X = [];
    for j = 1:numel(xb)-1
        [bl, p0] = Intersect(p1, p2, [xb(j),yb(j)], [xb(j+1), yb(j+1)]);
        if(bl)
            X = [X,p0(1)];
        end
    end
end

function [bubble, a] = VeloctyBiasX(bubble, a, x1, x2, factor)
    p = bubble.position; a = 0*a;
    if(x1 < p(1) && p(1) < x2)
        umax = norm(bubble.u);
        bubble.u = [factor*umax,0.5*bubble.u(2)];
    end
end

function [bubble, a] = VeloctyBiasY(bubble, a, y1, y2, factor)
end

function [bubble, a] = VeloctyBias2(bubble, a)
end