function [Wing] = GenerateSymMesh(Wing,Flap)

% Generation of the lattice around the wing. First, the panels are
% construct thanks to the discretization defined. Second, the control
% points are generated, where the aerodynamic forces take place. Finally,
% the normal vectors to the lattice surface.


%% DISCRETIZATION

b1 = Flap.Parameters.b1;
l1 = Flap.Parameters.l1;
b  = Wing.Parameters.b;
cf = Flap.Parameters.cf;
xc = 1 - cf/Wing.Geometry.Cr;

x0 = (Geom_dist(0,xc,Flap.Parameters.Nc+1,Wing.Parameters.bias_x,true))';
x1 = (Geom_dist(xc,1,Flap.Parameters.Nf+1,Wing.Parameters.bias_x,true))';

y1 = (Geom_dist(-b1,b1,...
    2*Flap.Parameters.Ns1+1,Wing.Parameters.bias_y,true))';
yf = (Geom_dist(b1,b1+l1,...
    Flap.Parameters.Nsf+1,Wing.Parameters.bias_y,true))';
y3 = (Geom_dist(b1+l1,b/2,...
    Flap.Parameters.Ns3+1,Wing.Parameters.bias_y,true))';

% Auxiliary surface index

S1 = Flap.Parameters.Ns3+1;
S2 = S1 + Flap.Parameters.Nsf;
S3 = S2 + 2*Flap.Parameters.Ns1;
S4 = S3 + Flap.Parameters.Nsf;
ny = S4 + Flap.Parameters.Ns3;
nx = Flap.Parameters.Nc + Flap.Parameters.Nf +1;
y  = zeros(ny,1);
for i = 1:ny

    if i <= S1
        y(i) = -y3(end-i+1);
    elseif i >= S1 && i <= S2
        y(i) = -yf(end-i+S1);
    elseif i > S2 && i <= S3
        y(i) = y1(i-S2+1);
    elseif i >= S3 && i <= S4
        y(i) = yf(i-S3+1);
    else
        y(i) = y3(i-S4+1);
    end
end

Wing.Mesh.Y = ones(nx,1)*y';
Wing.Mesh.X = zeros(nx,ny);
Wing.Mesh.Z = zeros(nx,ny);

Xs = zeros(nx,ny);
Zs = zeros(nx,ny);

for j = 1:ny
    for i = 1:Flap.Parameters.Nc+1
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

for j = 1:ny
    for i = Flap.Parameters.Nc+1:Flap.Parameters.Nc+Flap.Parameters.Nf+1
        yij = Wing.Mesh.Y(i,j);
        x1ij = x1(i - Flap.Parameters.Nc);
        c = real(Wing.Geometry.C(yij));
        xij = Wing.Geometry.Cr/4+(x1ij-0.25)*c+abs(yij)*...
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
        
        Xs(i,j) = abs(yij)*tan(Wing.Parameters.Sweepr)+[1,0]*(Ay*[xij - ...
            Wing.Geometry.Cr/4 - abs(yij)*tan(Wing.Parameters.Sweepr),zij]')...
            + Wing.Geometry.Cr/4;
        Zs(i,j) = [0,1]*(Ay*[xij-Wing.Geometry.Cr/4 - abs(yij)*...
            tan(Wing.Parameters.Sweepr),zij]')+abs(yij)*...
            tan(Wing.Parameters.dihedralr);
    end
end






%% CONTROL POINTS

Xc     = zeros(nx-1,ny-1);
Yc     = zeros(nx-1,ny-1);
Zc     = zeros(nx-1,ny-1);
dZc    = zeros(nx-1,ny-1);
deltaP = zeros(nx-1,ny-1);
Nvec   = zeros(nx-1,ny-1,3);

for j = 1:ny-1
    for i = 1:nx-1
        
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

        if i > Flap.Parameters.Nc && ((j >= S3 && j < S4) || (j >= S1 && j < S2))
            
            x1 = (Xc(i,j)-epsilon-tan(Wing.Parameters.Sweepr)*abs(Yc(i,j))...
            - Wing.Geometry.C(0)/4)/Wing.Geometry.C(Yc(i,j))+0.25;
        
            x2 = (Xc(i,j)+epsilon-tan(Wing.Parameters.Sweepr)*abs(Yc(i,j))-...
            Wing.Geometry.C(0)/4)/Wing.Geometry.C(Yc(i,j))+0.25;
            
            x1 = x1*cos(Flap.Parameters.dr);
            x2 = x2*cos(Flap.Parameters.dr);
            
            z1 = interp1([0,Wing.Parameters.b/2],...
            [Wing.Geometry.zc_r(x1),Wing.Geometry.zc_t(x1)],...
            abs(Yc(i,j)))*Wing.Geometry.C(Yc(i,j));
        
            z2 = interp1([0,Wing.Parameters.b/2],...
            [Wing.Geometry.zc_r(x2),Wing.Geometry.zc_t(x2)],...
            abs(Yc(i,j)))*Wing.Geometry.C(Yc(i,j));
        
            z1 = z1 - (Wing.Geometry.C(Yc(i,j)) - x1)*sin(Flap.Parameters.dr);
            z2 = z2 - (Wing.Geometry.C(Yc(i,j)) - x2)*sin(Flap.Parameters.dr);
            
            DeltaZ = z2-z1;
            dZc(i,j) = DeltaZ/DeltaX;
       
            deltaP(i,j) = atan(dZc(i,j));
            
            v1 = [Xc(i,j)*cos(Flap.Parameters.dr),Yc(i,j),Zc(i,j)-(Wing.Geometry.C(Yc(i,j)) - Xc(i,j))*sin(Flap.Parameters.dr)]...
                -[Xs(i,j)*cos(Flap.Parameters.dr),Wing.Mesh.Y(i,j),Zs(i,j)-(Wing.Geometry.C(Yc(i,j)) - Xc(i,j))*sin(Flap.Parameters.dr)];
            v2 = [Xc(i,j)*cos(Flap.Parameters.dr),Yc(i,j),Zc(i,j)-(Wing.Geometry.C(Yc(i,j)) - Xc(i,j))*sin(Flap.Parameters.dr)]...
                -[Xs(i,j+1)*cos(Flap.Parameters.dr),Wing.Mesh.Y(i,j+1),Zs(i,j+1)-(Wing.Geometry.C(Yc(i,j)) - Xc(i,j))*sin(Flap.Parameters.dr)];
        end
        
        nvec = cross(v2,v1)/norm(cross(v2,v1));
        
        A = [cos(deltaP(i,j)),0,-sin(deltaP(i,j));...
            0,1,0;...
            sin(deltaP(i,j)),0,cos(deltaP(i,j))];
        Nvec(i,j,:) = (A*nvec')';
    end
end

% Nvec(:,S1:S2,:) = fliplr(Nvec(:,S3:S4,:));



Xt = zeros(nx-1,ny);
Zt = zeros(nx-1,ny);

for i = 1:nx-1
    for j = 1:ny
        Xt(i,j) = 0.75*Xs(i,j)+0.25*Xs(i+1,j);
        Zt(i,j) = 0.75*Zs(i,j)+0.25*Zs(i+1,j);
    end
end



figure
hold on
mesh(Xs,Wing.Mesh.Y,Zs)
plot3(Xc,Yc,Zc,'.k')
quiver3(Xc,Yc,Zc,Nvec(:,:,1),Nvec(:,:,2),Nvec(:,:,3),0.5)
plot3(Xt',(Wing.Mesh.Y(1:nx-1,:))',Zt','--r')
hold off
axis('equal')
grid


Xt = [Xt; Xt(nx-1,:) + (Xt(nx-1,:) -  ...
    Xt(nx-1-1,:))*1000] ;
Zt = [Zt; Zt(nx-1,:) + (Zt(nx-1,:) - ...
    Zt(nx-1-1,:))*1000] ;

Wing.Mesh.Control.X = Xc;
Wing.Mesh.Control.Y = Yc;
Wing.Mesh.Control.Z = Zc;
Wing.Mesh.Node.X    = Xt;
Wing.Mesh.Node.Y    = Wing.Mesh.Y;
Wing.Mesh.Node.Z    = Zt;

Wing.Mesh.Xs        = Xs;
Wing.Mesh.Zs        = Zs;


Wing.Mesh.Nvec      = Nvec;

end %eof