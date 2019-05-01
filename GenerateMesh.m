function [Wing] = GenerateMesh(Wing)

% Generation of the lattice around the wing. First, the panels are
% construct thanks to the discretization defined. Second, the control
% points are generated, where the aerodynamic forces take place. Finally,
% the normal vectors to the lattice surface.

%% DISCRETIZATION


x0 = (Geom_dist(0,1,Wing.Parameters.nx,Wing.Parameters.bias_x,true))';

y = (Geom_dist(-Wing.Parameters.b/2,Wing.Parameters.b/2,...
    Wing.Parameters.ny,Wing.Parameters.bias_y,true))';
y(Wing.Parameters.Nss+1) = 0; 

Wing.Mesh.Y = ones(Wing.Parameters.nx,1)*y';
Wing.Mesh.X = zeros(Wing.Parameters.nx,Wing.Parameters.ny);
Wing.Mesh.Z = zeros(Wing.Parameters.nx,Wing.Parameters.ny);


for j = 1:Wing.Parameters.ny
    for i = 1:Wing.Parameters.nx
        yij = Wing.Mesh.Y(i,j);
        x0ij = x0(i);
        c = real(Wing.Geometry.C(yij));
        xij = Wing.Geometry.Cr/4+(x0ij-0.25)*c+abs(yij)*...
            tan(Wing.Parameters.Sweepr);
        zij = 0 ; 
        eps_y = Wing.Parameters.tor(yij);
        Ay = [cos(eps_y),sin(eps_y);-sin(eps_y),cos(eps_y)]; ...
            
        Wing.Mesh.X(i,j) = abs(yij)*tan(Wing.Parameters.Sweepr)+[1,0]*(Ay*...
            [xij-Wing.Geometry.Cr/4 - abs(yij)*...
            tan(Wing.Parameters.Sweepr),zij]')+ Wing.Geometry.Cr/4;
        Wing.Mesh.Z(i,j) = [0,1]*(Ay*[xij-Wing.Geometry.Cr/4 - abs(yij)*...
            tan(Wing.Parameters.Sweepr),zij]')+abs(yij)*...
            tan(Wing.Parameters.dihedral); 
    end
end



Xs = zeros(Wing.Parameters.nx,Wing.Parameters.ny);
Zs = zeros(Wing.Parameters.nx,Wing.Parameters.ny);
for j = 1:Wing.Parameters.ny
    for i = 1:Wing.Parameters.nx
        yij = Wing.Mesh.Y(i,j);
        x0ij = x0(i);
        c = real(Wing.Geometry.C(yij));
        xij = Wing.Geometry.Cr/4+(x0ij-0.25)*c+abs(yij)*...
            tan(Wing.Parameters.Sweepr);
        zij = 0;
        eps_y = Wing.Parameters.tor(yij);
        Ay = [cos(eps_y),sin(eps_y);-sin(eps_y),cos(eps_y)]; ...
            % matriz de giro para torsión
        Xs(i,j) = abs(yij)*tan(Wing.Parameters.Sweepr)+[1,0]*(Ay*[xij - ...
            Wing.Geometry.Cr/4 - abs(yij)*tan(Wing.Parameters.Sweepr),zij]')...
            + Wing.Geometry.Cr/4;
        Zs(i,j) = [0,1]*(Ay*[xij-Wing.Geometry.Cr/4 - abs(yij)*...
            tan(Wing.Parameters.Sweepr),zij]')+abs(yij)*...
            tan(Wing.Parameters.dihedralr);
        
        Wing.Geometry.torr(j) = eps_y;
        Wing.Geometry.tor(j)  = eps_y*180/pi;
    end
end


%% CONTROL POINTS

Xc     = zeros(Wing.Parameters.Nc,2*Wing.Parameters.Nss);
Yc     = zeros(Wing.Parameters.Nc,2*Wing.Parameters.Nss);
Zc     = zeros(Wing.Parameters.Nc,2*Wing.Parameters.Nss);
dZc    = zeros(Wing.Parameters.Nc,2*Wing.Parameters.Nss);
deltaP = zeros(Wing.Parameters.Nc,2*Wing.Parameters.Nss);
Nvec   = zeros(Wing.Parameters.Nc,2*Wing.Parameters.Nss,3);

for j = 1:2*Wing.Parameters.Nss
    for i = 1:Wing.Parameters.Nc
        
        Yc(i,j) = 0.5*(Wing.Mesh.Y(i,j)+Wing.Mesh.Y(i,j+1));
        Xc(i,j) = 0.5*(Xs(i,j)+Xs(i,j+1))+3/8*(Xs(i+1,j)-Xs(i,j)+...
            Xs(i+1,j+1)-Xs(i,j+1));
        Zc(i,j) = 0.5*(Zs(i,j)+Zs(i,j+1))+3/8*(Zs(i+1,j)-...
            Zs(i,j)+Zs(i+1,j+1)-Zs(i,j+1));
        
        epsilon = Wing.Geometry.C(Yc(i,j))*0.05; ...
            
        DeltaX = 2*epsilon;
        
        x1 = (Xc(i,j)-epsilon-tan(Wing.Parameters.Sweepr)*abs(Yc(i,j))...
            - Wing.Geometry.C(0)/4)/Wing.Geometry.C(Yc(i,j))+0.25;
        
        x2 = (Xc(i,j)+epsilon-tan(Wing.Parameters.Sweepr)*abs(Yc(i,j))-...
            Wing.Geometry.C(0)/4)/Wing.Geometry.C(Yc(i,j))+0.25;
        
        z1 = interp1([0,Wing.Parameters.b/2],...
            [Wing.Geometry.zc_r(x1),Wing.Geometry.zc_t(x1)],...
            abs(Yc(i,j)))*Wing.Geometry.C(Yc(i,j));
        
        z2 = interp1([0,Wing.Parameters.b/2],...
            [Wing.Geometry.zc_r(x2),Wing.Geometry.zc_t(x2)],...
            abs(Yc(i,j)))*Wing.Geometry.C(Yc(i,j));
        
        DeltaZ = z2-z1;
        dZc(i,j) = DeltaZ/DeltaX;
       
        deltaP(i,j) = atan(dZc(i,j));
        
        v1 = [Xc(i,j),Yc(i,j),Zc(i,j)]-[Xs(i,j),Wing.Mesh.Y(i,j),Zs(i,j)];
        v2 = [Xc(i,j),Yc(i,j),Zc(i,j)]-[Xs(i,j+1),Wing.Mesh.Y(i,j+1),Zs(i,j+1)];
        nvec = cross(v2,v1)/norm(cross(v2,v1));
        
        A = [cos(deltaP(i,j)),0,-sin(deltaP(i,j));...
            0,1,0;...
            sin(deltaP(i,j)),0,cos(deltaP(i,j))];
        Nvec(i,j,:) = (A*nvec')';
    end
end





Xt = zeros(Wing.Parameters.Nc,Wing.Parameters.ny);
Zt = zeros(Wing.Parameters.Nc,Wing.Parameters.ny);

for i = 1:Wing.Parameters.Nc
    for j = 1:Wing.Parameters.ny
        Xt(i,j) = 0.75*Xs(i,j)+0.25*Xs(i+1,j);
        Zt(i,j) = 0.75*Zs(i,j)+0.25*Zs(i+1,j);
    end
end



figure
hold on
mesh(Xs,Wing.Mesh.Y,Zs)
plot3(Xc,Yc,Zc,'.k')
quiver3(Xc,Yc,Zc,Nvec(:,:,1),Nvec(:,:,2),Nvec(:,:,3),0.5)
plot3(Xt',(Wing.Mesh.Y(1:Wing.Parameters.Nc,:))',Zt','--r')
hold off
axis('equal')
grid


Xt = [Xt; Xt(Wing.Parameters.Nc,:) + (Xt(Wing.Parameters.Nc,:) -  ...
    Xt(Wing.Parameters.Nc-1,:))*1000] ;
Zt = [Zt; Zt(Wing.Parameters.Nc,:) + (Zt(Wing.Parameters.Nc,:) - ...
    Zt(Wing.Parameters.Nc-1,:))*1000] ;


Wing.Mesh.X         = Wing.Mesh.X + Wing.Parameters.xoffset;
Wing.Mesh.Y         = Wing.Mesh.Y + Wing.Parameters.yoffset;
Wing.Mesh.Z         = Wing.Mesh.Z + Wing.Parameters.zoffset;
Wing.Mesh.Control.X = Xc + Wing.Parameters.xoffset;
Wing.Mesh.Control.Y = Yc + Wing.Parameters.yoffset;
Wing.Mesh.Control.Z = Zc + Wing.Parameters.zoffset;
Wing.Mesh.Node.X    = Xt + Wing.Parameters.xoffset;
Wing.Mesh.Node.Y    = Wing.Mesh.Y;
Wing.Mesh.Node.Z    = Zt + Wing.Parameters.zoffset;

Wing.Mesh.Xs        = Xs + Wing.Parameters.xoffset;
Wing.Mesh.Zs        = Zs + Wing.Parameters.zoffset;


Wing.Mesh.Nvec      = Nvec;