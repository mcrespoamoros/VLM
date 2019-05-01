function DisplayCase2(Wing,Tail)
 
C  = ones(size(Wing.Mesh.X));
C2 = ones(size(Tail.Mesh.X));
 
figure()
mesh(Wing.Mesh.X,Wing.Mesh.Y,Wing.Mesh.Z,C,'LineWidth',1); hold on
mesh(Tail.Mesh.X,Tail.Mesh.Y,Tail.Mesh.Z,C2,'LineWidth',1)
axis equal
colormap('vga')
xlabel('Chordwise [m]','Interpreter','latex','FontSize',15)
ylabel('Spanwise [m]','Interpreter','latex','FontSize',15)


