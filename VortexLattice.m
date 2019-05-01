function [Wing] = VortexLattice(Wing,FC)

% Defining control points

XC    = Wing.Mesh.Control.X;
YC    = Wing.Mesh.Control.Y;
ZC    = Wing.Mesh.Control.Z;

% Defining Nodes

Xn    = Wing.Mesh.Node.X;
Yn    = Wing.Mesh.Node.Y;
Zn    = Wing.Mesh.Node.Z;


% Auxiliar data

Xs    = Wing.Mesh.Xs;
Zs    = Wing.Mesh.Zs;
Nvec  = Wing.Mesh.Nvec;

% Definition of AIC

Aij = zeros(Wing.Parameters.Nt,Wing.Parameters.Nt);


for i1 = 1:Wing.Parameters.Nc
    for j1 = 1:2*Wing.Parameters.Nss
                    PC = [XC(i1,j1), YC(i1,j1), ZC(i1,j1)];
        for i2 = 1:Wing.Parameters.Nc
            for j2 = 1:2*Wing.Parameters.Nss
                
                    
                    P1 = [Xn(i2,j2), Yn(i2,j2), Zn(i2,j2)];
                    P2 = [Xn(i2,j2+1), Yn(i2,j2+1), Zn(i2,j2+1)];
                    P3 = [Xn(i2+1,j2+1), Yn(i2+1,j2+1), Zn(i2+1,j2+1)];
                    P4 = [Xn(i2+1,j2), Yn(i2+1,j2), Zn(i2+1,j2)];

                    
                    
                % Vortex Head 1 - 2
                
                r21 = P2 - P1;
                r1  = PC - P1;
                r2  = PC - P2;
                
           
                psi   = cross(r1,r2)/((norm(cross(r1,r2)))^2); 
                omega = r21*(r1/norm(r1)-r2/norm(r2))'; 
                V12   = 1/(4*pi)*psi*omega;
                
                if norm(cross(r1/norm(r1),r2/norm(r2))) < 1e-12
                    V12 = 0;
                end
                
                % Right side Vortex 2 - 3
                
                r32 = P3 - P2;
                r3 = PC - P3;

                
                psi = cross(r2,r3)/((norm(cross(r2,r3)))^2); 
                omega = r32*(r2/norm(r2)-r3/norm(r3))';
                V23 = 1/(4*pi)*psi*omega;
                
                 if norm(cross(r2/norm(r2),r3/norm(r3))) < 1e-12
                    V23 = 0;
                 end
                
                 
                % Free Stream Vortex 3 - 4
                
                
                r43 = P4 - P3;
                r4 = PC - P4;

                
                psi = cross(r3,r4)/((norm(cross(r3,r4)))^2); 
                omega = r43*(r3/norm(r3)-r4/norm(r4))'; 
                V34 = 1/(4*pi)*psi*omega;
                
                
                if norm(cross(r3/norm(r3),r4/norm(r4))) < 1e-12
                    V34 = 0;
                end
                
                
                % Left Side Vortex 4 - 1
                
                
                r14 = P1 - P4;

                psi = cross(r4,r1)/((norm(cross(r4,r1)))^2); 
                omega = r14*(r4/norm(r4)-r1/norm(r1))'; 
                V41 = 1/(4*pi)*psi*omega;
                if norm(cross(r4/norm(r4),r1/norm(r1))) < 1e-12
                    V41 = 0;
                end
                
                % Induced Velocity
                
                Vind = V12+V23+V34+V41;
                
                % Normal Component
                
                vn = Vind*[Nvec(i1,j1,1);Nvec(i1,j1,2);Nvec(i1,j1,3)];

                i = i1+Wing.Parameters.Nc*(j1-1);
                j = i2+Wing.Parameters.Nc*(j2-1);
                Aij(i,j) = Vind(end);


            end
        end
    end
end

Wing.VLM.AIC     = Aij;




for k = 1:FC.naoa
    
    % Struct name
    
    if FC.aoa(k) <0
        
        FC.aoa(k)         = -FC.aoa(k);
        name              = sprintf('AoAneg%1i', 100*FC.aoa(k));
        name              = strrep(name,'.','');
        FC.aoa(k)         = -FC.aoa(k);
    else
        
        name              = sprintf('AoA%1i', 100*FC.aoa(k));
        name              = strrep(name,'.','');
    end
        
    Wing.VLM.(name).clj   = zeros(1,2*Wing.Parameters.Nss);
    
    Talpha = [cos(FC.aoar(k)),0,-sin(FC.aoar(k));...
        0,1,0;...
        sin(FC.aoar(k)),0,cos(FC.aoar(k))];
    Tbeta  = [cos(FC.betar),sin(FC.betar),0;...
        -sin(FC.betar),cos(FC.betar),0;...
        0,0,1];
    
    U      = (FC.U*Tbeta*Talpha*[1,0,0]')';
    
    Wing.VLM.(name).Wyij   = zeros(Wing.Parameters.Nc*2*Wing.Parameters.Nss,1);
    
    for i1 = 1:Wing.Parameters.Nc
        for j1 = 1:2*Wing.Parameters.Nss
            vinf_n = U*[Nvec(i1,j1,1);Nvec(i1,j1,2);Nvec(i1,j1,3)];
            i = i1+Wing.Parameters.Nc*(j1-1);
            Wing.VLM.(name).Wyij(i) = -vinf_n;
        end
    end
    

    
    Wing.VLM.(name).Gamma = Aij\Wing.VLM.(name).Wyij;
    
    
    
    figures = true;
    
    
    %% Postprocessing
    

    
    Gamma2 = zeros(Wing.Parameters.Nc,2*Wing.Parameters.Nss);
    for j = 1:2*Wing.Parameters.Nss
        for i = 1:Wing.Parameters.Nc
            Gamma2(i,j) = Wing.VLM.(name).Gamma(i+Wing.Parameters.Nc*(j-1));
        end
    end
    
    Gamma2 = [Gamma2(1,:) ; Gamma2(2:end,:) - Gamma2(1:end-1,:)] ;
    
    if figures == true
        esc = 1/max(max(abs(Gamma2))); 
        figure
        % titlef = sprintf('\Gamma, aoa = %2.1f º',...
        %     FC.aoa(k));
        hold on
        surf(Xs,Wing.Mesh.Y,Zs,Gamma2)
        % colormap (jet)
        colorbar
        if Wing.Parameters.Nc == 1
            mesh ([Wing.Mesh.Control.X;Wing.Mesh.Control.X],...
                [Wing.Mesh.Control.Y;Wing.Mesh.Control.Y],...
                esc*[Gamma2;Gamma2],[Gamma2;Gamma2])
        else
            mesh (Wing.Mesh.Control.X,Wing.Mesh.Control.Y,esc*Gamma2,Gamma2,'facecol','none')
        end
        hold off
        az = 90;
        el = 90;
        view(az, el);
        % axis('equal')
        xlabel('Chordwise [m]')
        ylabel('Spanwise [m]')
        title('$\Gamma$')
        grid
        set(gca, 'DataAspectRatio', [1, 1, 1/Wing.Parameters.b])
    end
    
       
    Wing.VLM.(name).clij = zeros(Wing.Parameters.Nc,2*Wing.Parameters.Nss);
    Wing.VLM.(name).cdij = zeros(Wing.Parameters.Nc,2*Wing.Parameters.Nss);
    for i = 1:Wing.Parameters.Nc
        for j = 1:2*Wing.Parameters.Nss
            n = [Nvec(i,j,1),Nvec(i,j,2),Nvec(i,j,3)];
            cij = 0.5*(Xs(i+1,j)-Xs(i,j)+Xs(i+1,j+1)-Xs(i,j+1));
            Wing.VLM.(name).Cp = 2*Gamma2(i,j)/FC.U;
            Wing.VLM.(name).cdij(i,j) = Wing.VLM.(name).Cp/cij*U/FC.U*n';
%             Wing.VLM.(name).Cpij(i,j) = 2*Gamma2(i,j)/FC.U/cij;
            if Wing.VLM.(name).Cp == 0
                Wing.VLM.(name).clij(i,j) = 0;
            else
                Wing.VLM.(name).clij(i,j) = Wing.VLM.(name).Cp/...
                    abs(Wing.VLM.(name).Cp)*norm(Wing.VLM.(name).Cp/cij*n-...
                    Wing.VLM.(name).cdij(i,j)*U/FC.U);
            end
        end
    end
    
Wing.VLM.(name).Cpij = 2*Gamma2/FC.U;
   
    if figures == true
        esc2 = 1/max(max(abs(Wing.VLM.(name).clij)));
        figure
        
        hold on
        surf(Wing.Mesh.X,Wing.Mesh.Y,Wing.Mesh.Z,Wing.VLM.(name).clij)
        % colormap (jet)
        colorbar
        if Wing.Parameters.Nc == 1
            mesh ([Wing.Mesh.Control.X;Wing.Mesh.Control.X],...
                [Wing.Mesh.Control.Y;Wing.Mesh.Control.Y],...
                esc2*[Wing.VLM.(name).clij;Wing.VLM.(name).clij],...
                [Wing.VLM.(name).clij;Wing.VLM.(name).clij])
        else
            mesh(Wing.Mesh.Control.X,Wing.Mesh.Control.Y,...
                esc2*Wing.VLM.(name).clij,Wing.VLM.(name).clij,'facecol','none')
        end
        hold off
        az = 90;
        el = 90;
        view(az, el);
        % axis('equal')
        xlabel('Chordwise [m]')
        ylabel('Spanwise [m]')
        title('$C_{lij}$')
        grid
        set(gca, 'DataAspectRatio', [1, 1, 1/Wing.Parameters.b])
    end
    
    Z0     = zeros(Wing.Parameters.Nc,2*Wing.Parameters.Nss);
    
    if figures == true
        esc3 = 1/max(max(abs(Wing.VLM.(name).Cpij)));
        figure
%         titlef = sprintf('$C_{pij}$, AoA = %2.1f º',...
%             FC.aoa(k));
        hold on
        surf(Wing.Mesh.X,Wing.Mesh.Y,Wing.Mesh.Z,Wing.VLM.(name).Cpij)
        % colormap (jet)
        colorbar
        if Wing.Parameters.Nc == 1
            mesh ([Wing.Mesh.Control.X;Wing.Mesh.Control.X],...
                [Wing.Mesh.Control.Y;Wing.Mesh.Control.Y],...
                esc3*[Wing.VLM.(name).Cpij;Wing.VLM.(name).Cpij],...
                [Wing.VLM.(name).Cpij;Wing.VLM.(name).Cpij])
        else
            mesh(Wing.Mesh.Control.X,Wing.Mesh.Control.Y,...
                esc3*Wing.VLM.(name).Cpij,Wing.VLM.(name).Cpij,'facecol','none')
        end
        hold off
        az = 90;
        el = 90;
        view(az, el);
        % axis('equal')
        xlabel('Chordwise [m]')
        ylabel('Spanwise [m]')
        title('$C_{p}$')
        grid
        set(gca, 'DataAspectRatio', [1, 1, 1/Wing.Parameters.b])
    end
    
figure()
surf(Wing.Mesh.Control.X,Wing.Mesh.Control.Y,Gamma2); hold on
mesh(Wing.Mesh.Control.X,Wing.Mesh.Control.Y,Z0,'facecol','none')
colorbar
%colormap(jet)
az = 90;
el = 90;
view(az, el);
xlabel('Chordwise [m]')
ylabel('Spanwise [m]')
title('$\Gamma$')
grid on

figure()
surf(Wing.Mesh.Control.X,Wing.Mesh.Control.Y,Wing.VLM.(name).Cpij); hold on
mesh(Wing.Mesh.Control.X,Wing.Mesh.Control.Y,Z0,'facecol','none')
colorbar
%colormap(jet)
az = 90;
el = 90;
view(az, el);
xlabel('Chordwise [m]')
ylabel('Spanwise [m]')
title('$C_{p}$')
grid on

figure()
surf(Wing.Mesh.Control.X(:,Wing.Parameters.Nss+1:end),Wing.Mesh.Control.Y(:,Wing.Parameters.Nss+1:end),Wing.VLM.(name).Cpij(:,Wing.Parameters.Nss+1:end)); hold on
mesh(Wing.Mesh.Control.X(:,Wing.Parameters.Nss+1:end),Wing.Mesh.Control.Y(:,Wing.Parameters.Nss+1:end),Z0(:,Wing.Parameters.Nss+1:end),'facecol','none')
colorbar
%colormap(jet)
az = 90;
el = 90;
view(az, el);
xlabel('Chordwise [m]')
ylabel('Spanwise [m]')
title('$C_{p}$')
grid on    
    
    
    % Spanwise coefficient
    
    clij2 = zeros(Wing.Parameters.Nc,2*Wing.Parameters.Nss);
    cd    = zeros(Wing.Parameters.Nc,2*Wing.Parameters.Nss);
    
    for j = 1:2*Wing.Parameters.Nss
        for i = 1:Wing.Parameters.Nc
            cij = 0.5*(Xs(i+1,j)-Xs(i,j)+Xs(i+1,j+1)-Xs(i,j+1));
            clij2(i,j) = Wing.VLM.(name).clij(i,j)*cij/...
                Wing.Geometry.C(Wing.Mesh.Control.Y(i,j));
            cd(i,j) = Wing.VLM.(name).cdij(i,j)*cij/...
                Wing.Geometry.C(Wing.Mesh.Control.Y(i,j));
        end
        Wing.VLM.(name).clj(j) = sum(clij2(:,j));
        Wing.VLM.(name).cdy(j) = sum(cd(:,j));
    end

%% Plotting Cl

if figures == true
        figure
%         titlef = sprintf('$c_{l}$(y), Angle of Attack = %2.1f º', FC.aoa(k));
        plot([-Wing.Parameters.b/2,Wing.Mesh.Control.Y(1,:),Wing.Parameters.b/2],...
            [0,Wing.VLM.(name).clj,0])
        xlabel('Spanwise (m)')
        ylabel('$c_l(y)$')
        title('$C_{l}$')
        grid
end



figure()
plot(Wing.Mesh.Control.Y(1,Wing.Parameters.Nss+1:end),Wing.VLM.(name).clj(Wing.Parameters.Nss+1:end))
xlabel('Spanwise (m)')
ylabel('$c_l(y)$')
title('$C_{l}$')
grid on

%% Plotting Cp

c14   = floor(Wing.Parameters.Nc/4);
c34   = round(3*Wing.Parameters.Nc/4);
figure()
plot(Wing.Mesh.Control.Y(1,Wing.Parameters.Nss+1:end),Wing.VLM.(name).Cpij(1,Wing.Parameters.Nss+1:end)); hold on
plot(Wing.Mesh.Control.Y(c14,Wing.Parameters.Nss+1:end),Wing.VLM.(name).Cpij(c14,Wing.Parameters.Nss+1:end))
plot(Wing.Mesh.Control.Y(c34,Wing.Parameters.Nss+1:end),Wing.VLM.(name).Cpij(c34,Wing.Parameters.Nss+1:end));
xlabel('Spanwise (m)')
ylabel('$c_p$')
leg=legend('Leading Edge','Chord 1/4','Chord 3/4', 'Trailing Edge');
set(leg,'interpreter','latex','location','best')
grid on
title('$C_{p}$')

s14   = floor(Wing.Parameters.Nss+Wing.Parameters.Nss/4+1);
s34   = round(Wing.Parameters.Nss+3*Wing.Parameters.Nc/4+1);
figure()
plot(Wing.Mesh.Control.X(:,1),Wing.VLM.(name).Cpij(:,1)); hold on
plot(Wing.Mesh.Control.X(:,s14),Wing.VLM.(name).Cpij(:,s14))
plot(Wing.Mesh.Control.X(:,s34),Wing.VLM.(name).Cpij(:,s34));
xlabel('Chordwise (m)')
ylabel('$c_p$')
leg=legend('Root','1/4 Semi Span','3/4 Semi Span', 'Trailing Edge');
set(leg,'interpreter','latex','location','best')
grid on
title('$C_{p}$')




% Total Lift Coefficient
    
bj = zeros(1,2*Wing.Parameters.Nss);
bcj = zeros(1,2*Wing.Parameters.Nss);
for j = 1:2*Wing.Parameters.Nss
        bj(j) = Wing.Mesh.Y(1,j+1)-Wing.Mesh.Y(i,j);
        bcj(j) = bj(j)*Wing.Geometry.C(YC(1,j));
end
Wing.VLM.(name).cL  = 1/Wing.Parameters.Sw*Wing.VLM.(name).clj*bcj';
Wing.VLM.(name).cDi = 1/Wing.Parameters.Sw*Wing.VLM.(name).cdy*bcj';
Wing.VLM.(name).E   = Wing.VLM.(name).cL/Wing.VLM.(name).cDi;

% Result variables

CLplot(k)           = Wing.VLM.(name).cL;
CDplot(k)           = Wing.VLM.(name).cDi;
Eplot(k)            = Wing.VLM.(name).E;

%% Pitching moment
    
m0ij = zeros(Wing.Parameters.Nc,2*Wing.Parameters.Nss);
    
for j = 1:2*Wing.Parameters.Nss
   for i = 1:Wing.Parameters.Nc
            cij = 0.5*(Xs(i+1,j)-Xs(i,j)+Xs(i+1,j+1)-Xs(i,j+1));
            Xc4ij = 0.5*(Xn(i,j)+Xn(i,j+1));
            m0ij(i,j) = -Wing.VLM.(name).clij(i,j)*Xc4ij*cij*bj(j);
   end
end
Wing.VLM.(name).cM0y = 1/(Wing.Parameters.Sw*Wing.Geometry.C(0))*sum(sum(m0ij));


end

figure()
plot(Wing.Mesh.Y(1,Wing.Parameters.Nss+1:end),Wing.Geometry.tor(Wing.Parameters.Nss+1:end))
xlabel('Spanwise')
ylabel('Torsion (º)')
grid on
title('Torsion along Spanwise')

if FC.naoa > 1
    
    Cla  = CLplot./FC.aoar;
    
    figure()
    plot(FC.aoa,CLplot)
    xlabel('$\alpha$')
    ylabel('$C_L$')
    grid on
    title('$C_{L}$ vs $\alpha$')
    
    figure()
    plot(FC.aoa,CDplot)
    xlabel('$\alpha$')
    ylabel('$C_D$')
    grid on
    title('$C_{D}$ vs $\alpha$')
    
    figure()
    plot(FC.aoa,Eplot)
    xlabel('$\alpha$')
    ylabel('E')
    grid on
    title('Efficiency vs $\alpha$')
    
    figure()
    plot(CLplot,CDplot)
    xlabel('$C_L$')
    ylabel('$C_D$')
    grid on
    title('Polar')
    
end

end


