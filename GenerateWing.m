function [Wing]=GenerateWing(Wing)

if strcmp(Wing.Parameters.PF,'trap') == 1
    C = @(y) Ctrap(y,Wing);
elseif strcmp(Wing.Parameters.PF,'elps') == 1
    C = @(y) Celps(y,Wing);
end

if str2double(Wing.Parameters.prof_R) < 10000
    Wing.Geometry.zc_r = @(x) NACA4(Wing.Parameters.prof_R,x);
elseif (str2double(Wing.Parameters.prof_R) >= 10^4) && (str2double(Wing.Parameters.prof_R) < 10^5)
    Wing.Geometry.zc_r = @(x) NACA5(Wing.Parameters.prof_R,x);
else
    error('Introduce a valid NACA airfoil')
end

if str2double(Wing.Parameters.prof_T) < 10000
    Wing.Geometry.zc_t = @(x) NACA4(Wing.Parameters.prof_T,x);
elseif (str2double(Wing.Parameters.prof_T) >= 10^4) && (str2double(Wing.Parameters.prof_T) < 10^5)
    Wing.Geometry.zc_t = @(x) NACA5(Wing.Parameters.prof_T,x);
else
    error('Introduce a valid NACA airfoil')
end

% Wing chords

Wing.Geometry.C  = C;
Wing.Geometry.Cr = C(0);
Wing.Geometry.Ct = C(Wing.Parameters.b/2);


% Aerodynamic Center

Wing.Geometry.Xca = 0.25*C(0)+2*tan(Wing.Parameters.Sweepr)/...
    Wing.Parameters.Sw*simpson(@(y)C(y)*y,0,Wing.Parameters.b/2,10);
Wing.Geometry.Yca = 2/Wing.Parameters.Sw*simpson(@(y)C(y)*y,0,Wing.Parameters.b/2,10);



Wing.Geometry.CMG = Wing.Parameters.Sw/Wing.Parameters.b;
Wing.Geometry.CMA = 2/Wing.Parameters.Sw*simpson(@(y)(C(y))^2,0,Wing.Parameters.b/2,10);
