function S = simpson(f, li, ls, n)


x = linspace (li, ls, n);
S = 0;
for i = 1:n-1
    %puntos para la interpolacion parabolica
    x1 = x(i);
    x3 = x(i+1);
    x2 = (x1+x3)/2.;
    
    %Definicion de las matrices del sistema de ecs.
    X = [x1^2 x1 1; x2^2 x2 1; x3^2 x3 1];
    F = [f(x1); f(x2); f(x3)];
    %Solucion del sistema (coefs de la parabola)
    A = X\F;
    
    %Calculo del area
    I = A(1)/3.*(x3^3-x1^3)+A(2)/2.*(x3^2-x1^2)+A(3)*(x3-x1);
    %valor de la inteegral
    S = S + I;
end
%fprintf('Integral = %4.6f \n',S)


end

