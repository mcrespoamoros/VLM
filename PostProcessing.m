function PostProcessing(Wing,Tail,FC)

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
   
   
    %% Postprocessing
    figures = true;

    
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
        surf(Wing.Mesh.Xs,Wing.Mesh.Y,Wing.Mesh.Zs,Gamma2)
        surf(Tail.Mesh.Xs,Tail.Mesh.Y,Tail.Mesh.Zs,Gamma2)
        % colormap (jet)
        colorbar
        if Wing.Parameters.Nc == 1
            mesh ([Wing.Mesh.Control.X;Wing.Mesh.Control.X],...
                [Wing.Mesh.Control.Y;Wing.Mesh.Control.Y],...
                esc*[Gamma2;Gamma2],[Gamma2;Gamma2])
            mesh ([Tail.Mesh.Control.X;Tail.Mesh.Control.X],...
                [Tail.Mesh.Control.Y;Tail.Mesh.Control.Y],...
                esc*[Gamma2;Gamma2],[Gamma2;Gamma2])
        else
            mesh (Wing.Mesh.Control.X,Wing.Mesh.Control.Y,esc*Gamma2,Gamma2,'facecol','none')
            mesh (Tail.Mesh.Control.X,Tail.Mesh.Control.Y,esc*Gamma2,Gamma2,'facecol','none')
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
    
     
    if figures == true
        esc2 = 1/max(max(abs(Wing.VLM.(name).clij)));
        figure
        
        hold on
        surf(Wing.Mesh.X,Wing.Mesh.Y,Wing.Mesh.Z,Wing.VLM.(name).clij)
        surf(Tail.Mesh.X,Tail.Mesh.Y,Tail.Mesh.Z,Tail.VLM.(name).clij)
        % colormap (jet)
        colorbar
        if Wing.Parameters.Nc == 1
            mesh ([Wing.Mesh.Control.X;Wing.Mesh.Control.X],...
                [Wing.Mesh.Control.Y;Wing.Mesh.Control.Y],...
                esc2*[Wing.VLM.(name).clij;Wing.VLM.(name).clij],...
                [Wing.VLM.(name).clij;Wing.VLM.(name).clij])
            mesh ([Tail.Mesh.Control.X;Tail.Mesh.Control.X],...
                [Tail.Mesh.Control.Y;Tail.Mesh.Control.Y],...
                esc2*[Tail.VLM.(name).clij;Tail.VLM.(name).clij],...
                [Tail.VLM.(name).clij;Tail.VLM.(name).clij])
        else
            mesh(Wing.Mesh.Control.X,Wing.Mesh.Control.Y,...
                esc2*Wing.VLM.(name).clij,Wing.VLM.(name).clij,'facecol','none')
            mesh(Tail.Mesh.Control.X,Tail.Mesh.Control.Y,...
                esc2*Tail.VLM.(name).clij,Tail.VLM.(name).clij,'facecol','none')
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
        surf(Tail.Mesh.X,Tail.Mesh.Y,Tail.Mesh.Z,Tail.VLM.(name).Cpij)
        % colormap (jet)
        colorbar
        if Wing.Parameters.Nc == 1
            mesh ([Wing.Mesh.Control.X;Wing.Mesh.Control.X],...
                [Wing.Mesh.Control.Y;Wing.Mesh.Control.Y],...
                esc3*[Wing.VLM.(name).Cpij;Wing.VLM.(name).Cpij],...
                [Wing.VLM.(name).Cpij;Wing.VLM.(name).Cpij])
            mesh ([Tail.Mesh.Control.X;Tail.Mesh.Control.X],...
                [Tail.Mesh.Control.Y;Tail.Mesh.Control.Y],...
                esc3*[Tail.VLM.(name).Cpij;Tail.VLM.(name).Cpij],...
                [Tail.VLM.(name).Cpij;Tail.VLM.(name).Cpij])
        else
            mesh(Wing.Mesh.Control.X,Wing.Mesh.Control.Y,...
                esc3*Wing.VLM.(name).Cpij,Wing.VLM.(name).Cpij,'facecol','none')
            mesh(Tail.Mesh.Control.X,Tail.Mesh.Control.Y,...
                esc3*Tail.VLM.(name).Cpij,Tail.VLM.(name).Cpij,'facecol','none')
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
surf(Tail.Mesh.Control.X,Tail.Mesh.Control.Y,Gamma2)
mesh(Wing.Mesh.Control.X,Wing.Mesh.Control.Y,Z0,'facecol','none')
mesh(Tail.Mesh.Control.X,Tail.Mesh.Control.Y,Z0,'facecol','none')
colorbar
%colormap(jet)
az = 90;
el = 90;
view(az, el);
xlabel('Chordwise [m]')
ylabel('Spanwise [m]')
title('$\Gamma$')
grid on
hold off

figure()
surf(Wing.Mesh.Control.X,Wing.Mesh.Control.Y,Wing.VLM.(name).Cpij); hold on
surf(Tail.Mesh.Control.X,Tail.Mesh.Control.Y,Tail.VLM.(name).Cpij)
mesh(Wing.Mesh.Control.X,Wing.Mesh.Control.Y,Z0,'facecol','none')
mesh(Tail.Mesh.Control.X,Tail.Mesh.Control.Y,Z0,'facecol','none')
colorbar
%colormap(jet)
az = 90;
el = 90;
view(az, el);
xlabel('Chordwise [m]')
ylabel('Spanwise [m]')
title('$C_{p}$')
grid on
hold off

figure()
surf(Wing.Mesh.Control.X(:,Wing.Parameters.Nss+1:end),Wing.Mesh.Control.Y(:,Wing.Parameters.Nss+1:end),Wing.VLM.(name).Cpij(:,Wing.Parameters.Nss+1:end)); hold on
surf(Tail.Mesh.Control.X(:,Tail.Parameters.Nss+1:end),Tail.Mesh.Control.Y(:,Tail.Parameters.Nss+1:end),Tail.VLM.(name).Cpij(:,Tail.Parameters.Nss+1:end))
mesh(Wing.Mesh.Control.X(:,Wing.Parameters.Nss+1:end),Wing.Mesh.Control.Y(:,Wing.Parameters.Nss+1:end),Z0(:,Wing.Parameters.Nss+1:end),'facecol','none')
mesh(Tail.Mesh.Control.X(:,Tail.Parameters.Nss+1:end),Tail.Mesh.Control.Y(:,Tail.Parameters.Nss+1:end),Z0(:,Tail.Parameters.Nss+1:end),'facecol','none')
colorbar
%colormap(jet)
az = 90;
el = 90;
view(az, el);
xlabel('Chordwise [m]')
ylabel('Spanwise [m]')
title('$C_{p}$')
grid on    
    

end

end