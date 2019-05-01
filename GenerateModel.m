function [Model] = GenerateModel(Wing,Tail)

% This function joins the various aerodynamic surfaces in a single one

%% Model Parameters

Model.Parameters.Nc  = Wing.Parameters.Nc + Tail.Parameters.Nc;
Model.Parameters.Nss = Wing.Parameters.Nss;% + Tail.Parameters.Nss;
Model.Parameters.Nt  = Wing.Parameters.Nt + Tail.Parameters.Nt;

%% Mesh Parameters

Model.Mesh.X = [Wing.Mesh.X; Tail.Mesh.X];
Model.Mesh.Y = [Wing.Mesh.Y; Tail.Mesh.Y];
Model.Mesh.Z = [Wing.Mesh.Z; Tail.Mesh.Z];

% Control Points

Model.Mesh.Control.X = [Wing.Mesh.Control.X; Tail.Mesh.Control.X];
Model.Mesh.Control.Y = [Wing.Mesh.Control.Y; Tail.Mesh.Control.Y];
Model.Mesh.Control.Z = [Wing.Mesh.Control.Z; Tail.Mesh.Control.Z];

% Defining Nodes

Model.Mesh.Node.X    = [Wing.Mesh.Node.X; Tail.Mesh.Node.X];
Model.Mesh.Node.Y    = [Wing.Mesh.Node.Y; Tail.Mesh.Node.Y];
Model.Mesh.Node.Z    = [Wing.Mesh.Node.Z; Tail.Mesh.Node.Z];

% Auxiliar data

Model.Mesh.Xs        = [Wing.Mesh.Xs; Tail.Mesh.Xs];
Model.Mesh.Zs        = [Wing.Mesh.Zs; Tail.Mesh.Zs];
Model.Mesh.Nvec      = [Wing.Mesh.Nvec; Tail.Mesh.Nvec];

end
