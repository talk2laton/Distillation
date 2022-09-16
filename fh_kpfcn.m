function fh_kpfcn(H,E)
global HeaterTemp Heater Bubblef HeaterText Oil v3 v4
try
    switch E.Key
        case 'uparrow'
           HeaterTemp = min(2000, HeaterTemp+20);
        case 'downarrow'
           HeaterTemp = max(20, HeaterTemp-50);
    end
    Bubblef = (HeaterTemp-20)/1980;
    Heater.Color = Bubblef*[1,0,0];
    HeaterText.String = ['$',num2str(HeaterTemp),'^\circ C $'];
    Oil.FaceColor = [0,0,1] + Bubblef*0.5*([1,0,0]-[0,0,1]);
    v3 = repmat(Bubblef, size(v3)); v4 = repmat(Bubblef, size(v4));
    v3(end) = 1; v4(1) = 1;
    drawnow;
catch
end