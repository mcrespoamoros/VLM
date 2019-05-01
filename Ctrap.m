function C = Ctrap(y,Wing)
% Distribución de cuerdas pWing.Parameters.ARa ala trapecial

Wing.Parameters.b = (Wing.Parameters.AR*Wing.Parameters.Sw)^0.5;
Cr = 2/(Wing.Parameters.lambda+1)*Wing.Parameters.Sw/Wing.Parameters.b;
C = 2*Cr/Wing.Parameters.b*(Wing.Parameters.lambda-1)*abs(y)+Cr;

end