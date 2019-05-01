function C = Celps(y,Wing)


Wing.Parameters.b = (Wing.Parameters.AR*Wing.Parameters.Sw)^0.5;
Cr = 4*Wing.Parameters.Sw/(pi*b);
C = Cr*(1-(2*y/Wing.Parameters.b)^2)^0.5;

end

