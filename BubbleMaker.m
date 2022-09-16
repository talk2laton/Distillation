function Bubble = BubbleMaker(x,y,u,v,r,xs,ys)
atm = 1e5; gravity = 9.8; density = 1e3; gasconstant = 287;
Bubble.position = [x, y]; Bubble.u = [u, v]; 
Bubble.r = r; Bubble.collected = false; 
h = interp1(xs,ys,x);
p = atm + (h-y)*gravity*density;
Bubble.p = p; Bubble.m = p*4*pi*r^3/3/(gasconstant*300); 
t = 2*pi*linspace(0,1,21);
Bubble.bob = fill(x+r*cos(t), y+r*sin(t), 'w', EdgeAlpha = 0.5);
Bubble.Update = @(xs, ys, dt, components)update(xs, ys, dt, components);
    function bubble = update(xs, ys, dt, plthandles)
        h = interp1(xs,ys,Bubble.position(1));
        p = atm + max(0, h-Bubble.position(2))*gravity*density;
        r = Bubble.r*p/Bubble.p;
        F = 2e-4*gravity*density*4*pi*r^3/3;
        a = [0, F/Bubble.m]; 
        Bubble = movecollide(Bubble, a, dt, plthandles);
        Bubble.r = r; Bubble.p = p; 
        Bubble.bob.XData = Bubble.position(1)+r*cos(t);
        Bubble.bob.YData = Bubble.position(2)+r*sin(t);
        bubble = Bubble;
    end
end

function bubble = movecollide(bubble, a, dt, components)
    p1 = bubble.position;
    incomponent = 0;
    for i = 1:numel(components)
        component = components(i);
        X = component.Bound.XData; Y = component.Bound.YData;
        if(inpolygon(p1(1), p1(2), ...
                component.Fill.XData, component.Fill.YData))
            incomponent = 1;
            [bubble, a] = component.VelocityBias(bubble, a);
            v = bubble.u + a*dt;
            p2 = p1 + 0.5*(bubble.u + v)*dt;
            for j = 1:numel(X)-1
                [bl, p0] = Intersect(p1, p2, [X(j), Y(j)], [X(j+1), Y(j+1)]);
                if(bl)
                    dx = diff(X(j:j+1)); dy = diff(Y(j:j+1));
                    dl = hypot(dx, dy);
                    reflector = [dy,-dx]/dl;
                    dp = p2-p1; dp = dp/norm(dp);
                    c = dot(dp, reflector);
                    reflector = sign(c)*reflector;
                    dr = dp-2*reflector;
                    dr = dr/norm(dr);
                    v = norm(v)*dr;
                    p2 = p0 + norm(p2-p0)*dr;
                    break
                end
            end
            bubble.u = v;
            bubble.position = p2;
            return;
        end
        v = bubble.u + a*dt;
        p2 = p1 + 0.5*(bubble.u + v)*dt;
        for j = 1:numel(X)-1
            [bl, p0] = Intersect(p1, p2, [X(j), Y(j)], [X(j+1), Y(j+1)]);
            if(bl)
                dx = diff(X(j:j+1)); dy = diff(Y(j:j+1));
                dl = hypot(dx, dy);
                reflector = [dy,-dx]/dl;
                dp = p2-p1; dp = dp/norm(dp);
                c = dot(dp, reflector);
                reflector = sign(c)*reflector;
                dr = dp-2*reflector;
                dr = dr/norm(dr);
                v = norm(v)*dr;
                p2 = p0 + norm(p2-p0)*dr;
                break
            end
        end
        if(bl)
            break;
        end
    end
    bubble.collected = ~incomponent;
    bubble.u = v;
    bubble.position = p2;
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
