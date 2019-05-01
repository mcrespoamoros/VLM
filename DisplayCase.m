function DisplayCase(Wing,FC)


% This function generates the parameters needed to define the specified
% wing. For this, extreme points are defined in order to get the plant view
% of the wing. This points are going to be defined as the roots leading and
% trailing edge, and same for the tip leading and trailing edge. So for
% this, a series of points are generated with the subindex rl for root
% leading edge, rt root trailing, tl tip leading and tt tip trailing

for k=1:FC.naoa

    
    if FC.aoa(k) <0
        
        FC.aoa(k)         = -FC.aoa(k);
        name              = sprintf('AoAneg%1i', 100*FC.aoa(k));
        name              = strrep(name,'.','');
        FC.aoa(k)         = -FC.aoa(k);
    else
        
        name              = sprintf('AoA%1i', 100*FC.aoa(k));
        name              = strrep(name,'.','');
    end

Sweep    = Wing.Parameters.Sweepr;
Cr       = Wing.Parameters.Chord;
b        = Wing.Parameters.b/2;
lambda   = Wing.Parameters.lambda;
dihedral = Wing.Parameters.dihedralr;
aoa      = FC.aoar(k);

% Root leading edge

Wing.Geometry.(name).xrl = Wing.Parameters.xoffset;
Wing.Geometry.(name).yrl = Wing.Parameters.yoffset;
Wing.Geometry.(name).zrl = Wing.Parameters.zoffset;

% Tip leading edge

Wing.Geometry.(name).xtl = Wing.Geometry.(name).xrl + tan(Sweep)*b;
Wing.Geometry.(name).ytl = Wing.Geometry.(name).yrl + b;
Wing.Geometry.(name).ztl = Wing.Geometry.(name).zrl + tan(dihedral)*b;

% The leading edge of the wing is now defined

% Root trailing edge

Wing.Geometry.(name).xrt = Wing.Geometry.(name).xrl + Cr*cos(aoa);
Wing.Geometry.(name).yrt = Wing.Geometry.(name).yrl;
Wing.Geometry.(name).zrt = Wing.Geometry.(name).zrl - Cr*sin(aoa);

% Tip trailing edge

Wing.Geometry.(name).xtt = Wing.Geometry.(name).xtl + Cr*lambda*cos(aoa);
Wing.Geometry.(name).ytt = Wing.Geometry.(name).ytl;
Wing.Geometry.(name).ztt = Wing.Geometry.(name).zrl + tan(dihedral)*b - Cr*lambda*sin(aoa);

% In this function, a mesh is generated for evaluating the DLM. This mesh
% is saved inside the previously generated structure of Wing inside the
% Mesh substructure. 

% In this mesh, along with defining the Wing panels, the Control points are
% also defined and the Doublet lines. The control points are defined in
% their 3 coordinates as xC, yC and zC, with the index C referring as
% Control. 

% The Doublet Line is defined int the 1/4 chord of the panel. This line  
% will be defined as 3 limit points of this line. This points are
% defined by the limit points with the panel defined by the Albano & Rodden
% as the inner points, the outter points and the mid points. The notation
% used is xi for the x coordinate of the inner point, xo for the x
% coordinate of the outter point, and xm for the x coordinate of the mid
% point of the Doublet Linte. Idem for y and z.




% Creation of a 2D MESH (in x and y)

% For this mesh, the root chord and tip are divided in Nc+1 points, and
% same for the semi span coordinates, Ns+1, for the creation of the mesh
% nodes. There is always an N+1 node, where N is the number of panels. This
% nodes are defined equally spaced, so for this, knowing the limit
% coordinates defined in Wing.Geometry, a linspace is used to quickly
% generate this nodes coordinates. 

Nc      = Wing.Parameters.Nc;
Ns      = Wing.Parameters.Nss;
Nt      = Wing.Parameters.Nt;
torsion = Wing.Parameters.theta0r;
Root_x  = linspace(Wing.Geometry.(name).xrl,Wing.Geometry.(name).xrt,Nc+1); 
Tip_x   = linspace(Wing.Geometry.(name).xtl,Wing.Geometry.(name).xtt,Nc+1);
Span    = linspace(Wing.Geometry.(name).yrl,Wing.Geometry.(name).ytl,Ns+1);
Root_z  = linspace(Wing.Geometry.(name).zrl,Wing.Geometry.(name).zrt,Nc+1);
Tip_z   = linspace(Wing.Geometry.(name).ztl,Wing.Geometry.(name).ztt,Nc+1);

% Now, the full X, Y and Z coordinates are defined of the whole Wing. mx
% and mz are the slope used for the calculation of the points as a function
% (y-y0)=m(x-x0). Knowing the x0, y and y0 points, the x is defined.

for i=1:(Nc+1)
   
   mx= (Root_x(i)-Tip_x(i))/(Span(1)-Span(end));
   mz= (Root_z(i)-Tip_z(i))/(Span(1)-Span(end));
   
   for j=1:(Ns+1)
      
       X(i,j) = -(Span(1)-Span(j))*mx+Root_x(i);
       Y(i,j) =  Span(j);
       Z(i,j) = -(Span(1)-Span(j))*mz+Root_z(i);
       tor    = torsion*(1-Span(1,j)^2/(2*Span(end)^2));
       xo     = Cr/2 - X(i,j);
       Z(i,j) = Z(i,j) + sin(tor)*xo;
      
   end
    
end


for ii=1:(Nc+1)
    
    for jj=1:Ns
        
    Xmid(ii,jj) = (X(ii,jj+1)+X(ii,jj))/2;
    Ymid(ii,jj) = (Y(ii,jj+1)+Y(ii,jj))/2;
    Zmid(ii,jj) = (Z(ii,jj+1)+Z(ii,jj))/2;
    
    end
    
end

for i=1:Nc
    
    for j=1:Ns
        
        %INBOARD COORDINATES
        chord_i   = X(i+1,j)-X(i,j);
        xi(i,j)   = X(i,j)+chord_i/4;
        yi(i,j)   = Y(i,j);
        zi(i,j)   = Z(i,j);
        
        %OUTBOARD COORDINATES
        chord_o   = X(i+1,j+1)-X(i,j+1);
        xo(i,j)   = X(i,j+1)+chord_o/4;
        yo(i,j)   = Y(i,j+1);
        zo(i,j)   = Z(i,j+1);
        
        %MIDSPAN COORDINATES
        chord_m   = Xmid(i+1,j)-Xmid(i,j);
        xm(i,j)   = Xmid(i,j)+chord_m/4;
        ym(i,j)   = Ymid(i,j);
        zm(i,j)   = Zmid(i,j);
        
        %CONTROL POINTS
        xC(i,j)     = Xmid(i,j)+(3/4)*chord_m;
        yC(i,j)     = Ymid(i,j);
        zC(i,j)     = Zmid(i,j);
                

    end
    
end

xi        = reshape(xi',Nc*Ns,1);
yi        = reshape(yi',Nc*Ns,1);
zi        = reshape(zi',Nc*Ns,1);
xo        = reshape(xo',Nc*Ns,1);
yo        = reshape(yo',Nc*Ns,1);
zo        = reshape(zo',Nc*Ns,1);
xm        = reshape(xm',Nc*Ns,1);
xmid      = xm;
ym        = reshape(ym',Nc*Ns,1);
zm        = reshape(zm',Nc*Ns,1);




C  = ones(size(X));



figure()
mesh(X,Y,Z,C,'LineWidth',1); hold on
mesh(X,-Y,Z,C,'LineWidth',1)
axis equal
colormap('vga')
xlabel('Chordwise ','Interpreter','latex','FontSize',18)
ylabel('Spanwise','Interpreter','latex','FontSize',18)


figure()
plot3(xC(1,:),yC(1,:),zC(1,:),'bo','LineWidth',1); hold on
for i=1:Nt/2
   dx=[xo(i) xi(i)];  
   dy=[yo(i) yi(i)]; 
   dz=[zo(i) zi(i)]; 
   plot3(dx,dy,dz,'-r','LineWidth',1)
   dx=[xo(i) xi(i)];  
   dy=[-yo(i) -yi(i)]; 
   dz=[zo(i) zi(i)]; 
   plot3(dx,dy,dz,'-r','LineWidth',1)
end
mesh(X,Y,Z,C,'LineWidth',1);
plot3(xC,yC,zC,'bo','LineWidth',1) 
mesh(X,-Y,Z,C,'LineWidth',1)
plot3(xC,-yC,zC,'bo','LineWidth',1)

colormap('vga')
axis equal
set(gca, 'YAxisLocation', 'right')
xlabel('Chordwise [m]','Interpreter','latex','FontSize',18)
ylabel('Spanwise [m]','Interpreter','latex','FontSize',18)
% legend('Puntos de Control','Cabeza de torbellinos');

grid on


% i = round(Nc/2);
% dx=[X(i,i) X(i,i+1)];  
% dy=[Y(i,i) Y(i,i+1)]; 
% dz=[Z(i,i) Z(i,i+1)];
% dxi=[X(i,i) X(Nc+1,i)];
% dyi=[Y(i,i) Y(Nc+1,i)];
% dzi=[Z(i,i) Z(Nc+1,i)];
% dxo=[X(i,i) X(Nc,i+1)];
% dyo=[Y(i,i) Y(Nc,i+1)];
% dzo=[Z(i,i) Z(Nc,i+1)];
% dxfi=[X(Nc+1,i)   1.5*X(Nc+1,i)];
% dyfi=[Y(Nc+1,i)   Y(Nc+1,i)];
% dzfi=[Z(Nc+1,i)   Z(Nc+1,i)];
% dxfo=[X(Nc+1,i+1) 1.5*X(Nc+1,i+1)];
% dyfo=[Y(Nc+1,i+1) Y(Nc+1,i+1)];
% dzfo=[Z(Nc+1,i+1) Z(Nc+1,i+1)];
% figure()
% plot3(xC(i,i),yC(i,i),zC(i,i),'bo','LineWidth',1); hold on
% plot3(dx,dy,dz,'-r','LineWidth',1)
% plot3(xC,yC,zC,'bo','LineWidth',1); hold on
% plot3(xC,-yC,zC,'bo','LineWidth',1);
% mesh(X,Y,Z,C,'LineWidth',1);
% mesh(X,-Y,Z,C,'LineWidth',1);
% plot3(dxi,dyi,dzi,'^-r','LineWidth',1)
% plot3(dxo,dyo,dzo,'^-r','LineWidth',1)
% plot3(dxfi,dyfi,dzfi,'-r','LineWidth',1)
% plot3(dxfo,dyfo,dzfo,'-r','LineWidth',1)
% plot3(dxi,-dyi,dzi,'^-r','LineWidth',1)
% plot3(dxo,-dyo,dzo,'^-r','LineWidth',1)
% plot3(dxfi,-dyfi,dzfi,'-r','LineWidth',1)
% plot3(dxfo,-dyfo,dzfo,'-r','LineWidth',1)
% 
% colormap('vga')
% axis equal
% set(gca, 'YAxisLocation', 'right')
% xlabel('Chordwise [m]','Interpreter','latex','FontSize',15)
% ylabel('Spanwise [m]','Interpreter','latex','FontSize',15)
% legend('Control Points','Horseshoe Vortex');
% grid on

clear Xmid Ymid Zmid xi yi zi xo yo zo xm ym zm xmid
end

end
