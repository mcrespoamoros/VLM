function [zc, ze, zi] = NACA5(NACA, x)
%Puntos para perfiles NACA de 5 cifras
%   NACA = 'XXXXX'

A = str2double(NACA(1));
B = str2double(NACA(2));
m = B*0.05;
C = str2double(NACA(3));
t = str2double(NACA(4:5))/100;

% Línea media

if C == 0 % Línea media con curvatura simple
    if B == 1
        r = 0.058;
        k1 = A*361.4/2;
    elseif B == 2
        r = 0.126;
        k1 = A*51.64/2;
    elseif B == 3
        r = 0.2025;
        k1 = A*15.957/2;
    elseif B == 4
        r = 0.29;
        k1 = A*6.643/2;
    elseif B == 5
        r = 0.391;
        k1 = A*3.23/2;
    else
        error('El segundo dígito debe estar esntre 1 y 5, ambos incluidos')
    end
    if x < r
        zc = k1/6*(x^3-3*r*x^2+r^2*(3-r)*x);
    else
        zc = k1*r^3/6*(1-x);
    end
elseif C == 1 % Línea media "reflex" (doble curvatura)
    if B == 2
        r = 0.13;
        k1 = A*51.99/2;
        k2 = k1*0.000764;
    elseif B == 3
        r = 0.217;
        k1 = A*15.793/2;
        k2 = k1*0.00677;
    elseif B == 4
        r = 0.318;
        k1 = A*6.520/2;
        k2 = k1*0.0303;
    elseif B == 5
        r = 0.441;
        k1 = A*3.191/2;
        k2 = k1*0.1355;
    else
        error('El segundo dígito debe estar esntre 2 y 5, ambos incluidos')
    end
    if x < r
        zc = k1/6*((x-r)^3-k2/k1*(1-r)^3*x-r^3*x+r^3);
    else
        zc = k1/6*(k2/k1*(x-r)^3-k2/k1*(1-r)^3*x-r^3*x+r^3);
    end
else
    error('El tercer dígito debe ser 0 o 1')
end

%Espesor

a = 0.2969;
b = -0.1260;
c = -0.3516;
d = 0.2843;
e = -0.1015;

zt = 5*t*(a*sqrt(x)+b*x+c*x^2+d*x^3+e*x^4);

%Extrados

ze = zc+zt;

%intrados

zi = zc-zt;

end

