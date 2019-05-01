function y = Geom_dist(x0, x1, n, bias, sym)

    % x0: punto inicial
    % x1; punto final, > x0
    % n: número de nodos
    % bias: relacion entre el segmento fianl y el inicial (o central si
    %   es simétrica), por defecto = 1 (lineal)
    % sym: distribución simétrica (true), o normal (false), por defecto
    %   false
    
    % Valores por defecto
    if not(exist('bias','var'))
        bias = 1;
    end
    if not(exist('sym','var'))
        sym = false;
    end
    
    % Cálculo de la razón geométrica
    if sym
        r = bias^(1/(floor(n/2 - 1)));
    else
        r = bias^(1/(n-2));
    end
    
    if r == 1
        % distribución lineal
        Y = linspace(0,1,n);
    elseif not(sym)
        % homogeneamente creciente o decreciente
        Y = zeros(1,n);
        d = (1-r)/(1-r^(n-1));
        for i = 2:n
            Y(i) = Y(i-1)+d*r^(i-2);
        end
    else
        % simétrica
        Y = zeros(1,n);
        num = 1 - r;
        den = 2 - 2*r^((n-mod(n,2))/2) - mod(n+1,2)*(1 - r);
        d = num/den;
        for i = 2:n;
            m = floor(abs(i - 1 - n/2));
            Y(i) = Y(i-1) + d*r^m;
        end
    end
    y = Y.*(x1-x0)+x0;
end