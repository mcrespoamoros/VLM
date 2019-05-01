%% VORTEX LATTICE METHOD

% VLM for a wing with sweep, taper, angle of attack and torsion
% The sweep can be modeled for a non sweep wing (straigth wing), bigger 
% than 0 (conventional sweep) or smaller than 0 (regresive sweep)

% Additionally, a control element can be set as a Flap or Aileron. Being
% the symmetric case the Flap and the antisymmetric case the Aileron

% Crespo Amorós, Miguel - May 2018
% Last version - April 2019

%%

clear all; close all;

setPlot;
%% MAIN PARAMETERS

% Wing parameter definition

Wing.Parameters.Chord     = 5;
Wing.Parameters.AR        = 5;
Wing.Parameters.Sweep     = 45;
Wing.Parameters.Sweepr    = deg2rad(Wing.Parameters.Sweep);
Wing.Parameters.Sw        = Wing.Parameters.AR*Wing.Parameters.Chord^2;
Wing.Parameters.lambda    = 1;
Wing.Parameters.b         = (Wing.Parameters.AR*Wing.Parameters.Sw)^0.5;
Wing.Parameters.dihedral  = 0;
Wing.Parameters.dihedralr = deg2rad(Wing.Parameters.dihedral);

Wing.Parameters.xoffset   = 0;            
Wing.Parameters.yoffset   = 0;            
Wing.Parameters.zoffset   = 0;            

Wing.Parameters.prof_R  = '0000';
Wing.Parameters.prof_T  = '0000';
Wing.Parameters.PF      = 'trap';
Wing.Parameters.tor_d   = 0;%-0.25;
Wing.Parameters.theta0  = 0;
Wing.Parameters.theta0r = deg2rad(Wing.Parameters.theta0);
Wing.Parameters.tor     = @(y) Wing.Parameters.theta0r*(1-y^2/Wing.Parameters.b^2);

Wing.Parameters.Nc      = 10;
Wing.Parameters.nx      = Wing.Parameters.Nc+1;
Wing.Parameters.Nss     = 20;
Wing.Parameters.ny      = Wing.Parameters.Nss*2+1;
Wing.Parameters.Nt      = Wing.Parameters.Nc*Wing.Parameters.Nss*2;
Wing.Parameters.bias_x  = 1;
Wing.Parameters.bias_y  = 1;

% Tail parameter generation

Tail.Parameters.Chord     = 3;
Tail.Parameters.AR        = 6;
Tail.Parameters.Sweep     = 20;
Tail.Parameters.Sweepr    = deg2rad(Tail.Parameters.Sweep);
Tail.Parameters.Sw        = Tail.Parameters.AR*Tail.Parameters.Chord^2;
Tail.Parameters.lambda    = 0.8;
Tail.Parameters.b         = (Tail.Parameters.AR*Tail.Parameters.Sw)^0.5;
Tail.Parameters.dihedral  = 10;
Tail.Parameters.dihedralr = deg2rad(Tail.Parameters.dihedral);

Tail.Parameters.xoffset   = 20;            
Tail.Parameters.yoffset   = 0;            
Tail.Parameters.zoffset   = 5;            

Tail.Parameters.prof_R  = '0000';
Tail.Parameters.prof_T  = '0000';
Tail.Parameters.PF      = 'trap';
Tail.Parameters.tor_d   = -0.25;
Tail.Parameters.theta0  = 1;
Tail.Parameters.theta0r = deg2rad(Tail.Parameters.theta0);
Tail.Parameters.tor     = @(y) Tail.Parameters.theta0r*(1-y^2/Tail.Parameters.b^2);

Tail.Parameters.Nc      = 10;
Tail.Parameters.nx      = Tail.Parameters.Nc+1;
Tail.Parameters.Nss     = 10;
Tail.Parameters.ny      = Tail.Parameters.Nss*2+1;
Tail.Parameters.Nt      = Tail.Parameters.Nc*Tail.Parameters.Nss*2;
Tail.Parameters.bias_x  = 1;
Tail.Parameters.bias_y  = 1;

% Flap Parameters generation

Flap.Parameters.b1 = 5;
Flap.Parameters.cf = 1;
Flap.Parameters.l1 = 5;
Flap.Parameters.d  = -5;
Flap.Parameters.dr = deg2rad(Flap.Parameters.d);


Flap.Parameters.Ns1 = 5;
Flap.Parameters.Nsf = 10;
Flap.Parameters.Ns3 = 5;

Flap.Parameters.Nc  = 5;
Flap.Parameters.Nf  = 5;

% Flight Conditions

FC.aoa     = [2.1, 4.2, 6.3, 8.4, 10.5];
FC.aoar    = deg2rad(FC.aoa);
FC.beta    = 0;
FC.betar   = deg2rad(FC.beta);
FC.M       = 0.3;
FC.H       = 0;
[T,a,P,rho]= atmosisa(FC.H);
FC.U       = a*FC.M;
FC.q       = 0.5*rho*FC.U^2;
FC.naoa    = size(FC.aoar,2);

%% Wing Generation

casetype   = 'Standard';

[Wing]     = GenerateWing(Wing);

switch casetype
    
    case 'Standard'
        
        
        [Wing]     = GenerateMesh(Wing);
        DisplayCase(Wing,FC);
        [Wing]     = VortexLattice(Wing,FC);

    case 'Symmetric'

        [Wing]     = GenerateSymMesh(Wing,Flap);
        surf(Wing.Mesh.X,Wing.Mesh.Y,Wing.Mesh.Z);
        axis equal
        Wing.Parameters.Nc  = Flap.Parameters.Nc+Flap.Parameters.Nf;
        Wing.Parameters.Nss = Flap.Parameters.Ns1+Flap.Parameters.Nsf+Flap.Parameters.Ns3;
        Wing.Parameters.Nt  = Wing.Parameters.Nc*Wing.Parameters.Nss*2;

        [Wing]     = VortexLattice(Wing,FC);
        
    case 'AntiSymmetric'

        [Wing]     = GenerateAntMesh(Wing,Flap);
        surf(Wing.Mesh.X,Wing.Mesh.Y,Wing.Mesh.Z);
        axis equal
        Wing.Parameters.Nc  = Flap.Parameters.Nc+Flap.Parameters.Nf;
        Wing.Parameters.Nss = Flap.Parameters.Ns1+Flap.Parameters.Nsf+Flap.Parameters.Ns3;
        Wing.Parameters.Nt  = Wing.Parameters.Nc*Wing.Parameters.Nss*2;

        [Wing]     = VortexLattice(Wing,FC);
            
end

% [Wing]     = GenerateWing(Wing);
% [Wing]     = GenerateMesh(Wing);
% 
% % [Tail]     = GenerateWing(Tail);
% % [Tail]     = GenerateMesh(Tail);
% 
% DisplayCase(Wing,FC);
% % DisplayCase(Tail,FC);
% % DisplayCase2(Wing,Tail);
% 
% %% Full Model Generation
% 
% % [Model]    = GenerateModel(Wing,Tail);
% 
% %% Vortex Lattice Simulation
% 
% [Wing]     = VortexLattice(Wing,FC);
% % [Tail]     = VortexLattice(Tail,FC);
% 
% % PostProcessing(Wing,Tail,FC);
