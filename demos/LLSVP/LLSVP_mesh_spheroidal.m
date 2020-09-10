%% Build an Earth model containing an LLSVP.
% To be run from PlanetaryModels/demos/LLSVP.

%% Set up: Define parameters.

% Clear work space.
clear all;
clc;

% Choose whether to have an LLSVP (1) or not (0).
%llsvp = 0;

% Load input parameters.
fid = fopen('input_spheroidal.txt', 'r');
tet_max_vol = fscanf(fid, '%e', 1);
pOrder = fscanf(fid, '%d', 1);
fclose(fid);

% Define paths to modules.
addpath('../../modelbuilder/');
addpath('../../make_unit_spheres');
tetgen = '../../packages/tetgen1.5.0/tetgen';
distmesh = '../../packages/distmesh';

% Load radial Earth model.
% RD        Radii of discontinuity and indices of discontinuities.
% nlayer    Number of layers between discontinuities (3).
% MI        Model information: radius, density, Vp and Vs.
load ../../radialmodels/prem3L_noocean.mat

% Load radii of discontinuity.
R1 = RD(1,1);
R2 = RD(2,1);
R3 = RD(3,1);

% Exaggerate the ellipticity (set to NaN if not).
% Useful for plotting.
exagg = NaN;

% Load the ellipticity profile.
fid = fopen('ellipticity_profile.txt', 'r');
ellipticity_profile = fscanf(fid, '%e', [2,Inf]);
fclose(fid);
% Exaggerate ellipticity.
if ~isnan(exagg)
    ellipticity_profile(2, :) = ellipticity_profile(2, :)*exagg;
end
%
ellipticity_profile(1, :) = ellipticity_profile(1, :)*1.0E-3; % m to km.
r_ellipticity = ellipticity_profile(1, :);
ellipticity = ellipticity_profile(2, :);
%
% Remove duplicate points (doesn't change profile as ellipticity is continuous).
% The duplicate points break the interpolation (it works in Python...).
[r_ellipticity, i_unique, ~] = unique(r_ellipticity);
ellipticity = ellipticity(i_unique);
ellipticity_profile = ellipticity_profile(:, i_unique);

% Find the ellipticity values at the discontinuity radii.
eps_1 = interp1(r_ellipticity, ellipticity, R1);
eps_2 = interp1(r_ellipticity, ellipticity, R2);
eps_3 = interp1(r_ellipticity, ellipticity, R3);

% Choose finite element order (choose 1 or 2)
%pOrder  = 2;

% Set the maximum volume of tetrahedral element (km3).
%tet_max_vol = 1.0E7; % 439.4.
%tet_max_vol = 5.0E6; % 348.8.
%tet_max_vol = 1.0E6; % 204.0.

l_max = (2.0^(1.0/2.0))*(3.0^(1.0/3.0))*(tet_max_vol^(1.0/3.0)); 

% Set the ratio between the edge length of the tetrahedra at the free
% surface and at the middle of the mantle. The target edge length will
% increase linearly from the surface to the middle, then increase
% linearly again to the base of the mantle.
edge_length_ratio_mid_mantle_to_surf = 2.0;
tet_max_vol_surf = tet_max_vol/(edge_length_ratio_mid_mantle_to_surf^3.0);

% Calculate the number of triangles on the free surface.
a_surf          = 4.0*pi*(R1^2.0);
a_face_tet_max  = 8.0*(3.0^(-7.0/6.0))*(tet_max_vol_surf^(2.0/3.0));
n_tri_surf      = round(a_surf/a_face_tet_max);

% Calculate the number of triangles on the inner interfaces.
n_tri_cmb = round(n_tri_surf*((R2/R1)^2.0));
n_tri_icb = round(n_tri_surf*((R3/R1)^2.0));

%% Create or load triangulations on spheroidal interfaces.

tic;

% Calculate the edge-length parameter (d) for the distmesh function. This
% parameter controls the number of triangles in the surface meshing:
% approximately, n = k/(d**2) where k = 37.67795182. This approximation is
% accurate to within 10% for small meshes, and the accuracy improves for
% larger ones.
sqrtk   = 6.138236866;
d_surf  = sqrtk/sqrt(n_tri_surf);
d_cmb   = sqrtk/sqrt(n_tri_cmb);
d_icb   = sqrtk/sqrt(n_tri_icb);
%
% A bug in distmeshsurface can cause it to fail if the initial grid is not
% an exact multiple of the d parameter (initial spacing). Therefore, 
% we change the d parameter so it is equal to the nearest integer divisor.
d_grid_width    = 2.2;
d_surf          = d_grid_width/(round(d_grid_width/d_surf));
d_cmb           = d_grid_width/(round(d_grid_width/d_cmb));
d_icb           = d_grid_width/(round(d_grid_width/d_icb));

% Calculate unit spheroids (if they do not already exist)
dir_unit_spheroid = '../../output/unit_spheroids/';
unit_spheroid_file_surf   = unit_spheroid(d_surf, eps_1,   dir_unit_spheroid, distmesh);
unit_spheroid_file_cmb    = unit_spheroid(d_cmb,  eps_2,   dir_unit_spheroid, distmesh);
unit_spheroid_file_icb    = unit_spheroid(d_icb,  eps_3,   dir_unit_spheroid, distmesh);

% Load unit spheroids and then scale by radius.
% The 'unit spheroids' have semi-major axis of 1.0.
% The .mat files in ../../unitspheres/data list the number of triangles
% on the surface of the sphere.
% Can visualise with e.g. trimesh(t1, p1(:, 1), p1(:, 2), p1(:, 3)).
% p1    3D coordinates of triangulation points.
% np1   Number of triangulation points.
% t1    Triangulation indices.
% nt1   Number of triangles.
% Similarly for p2, np2, t2, nt2 and p3, np3, t3, nt3.

% Surface of the Earth.
load(unit_spheroid_file_surf)
R1_equator = R1*(1.0 + (eps_1/3.0));
p1              = R1_equator*p;         % Points on surface of Earth.
np1             = size(p1,1);
t1              = t;
nt1             = size(t1,1);

% Core-mantle boundary.
load(unit_spheroid_file_cmb)
R2_equator = R2*(1.0 + (eps_1/3.0));
p2              = R2_equator*p;         % Points on surface of Earth.
np2             = size(p2,1);
t2              = t;          
nt2             = size(t2,1);

% Surface of the inner core.
load(unit_spheroid_file_icb)
R3_equator = R3*(1.0 + (eps_1/3.0));
p3              = R3_equator*p;         % Points on surface of Earth.
np3             = size(p3,1);
t3              = t;
nt3             = size(t3,1);

% Merge points and triangulations.
% First, shift indices in triangulation.
t2      = t2 + np1;
t3      = t3 + np1 + np2;
%
pin     = cat(1, p1, p2, p3);
tin     = cat(1, t1, t2, t3);

% Define the name of the output directory and model file.
mod_name = 'mod';
dir_name = sprintf('spheroid_%06.1f_%4.2f_%1d',   ...
                            l_max, edge_length_ratio_mid_mantle_to_surf, pOrder);            
%dir_name_llsvp = sprintf('llsvp_%s', dir_name);
if ~isnan(exagg)
    dir_name_llsvp = sprintf('llsvp_exagg_%s', dir_name);
else
    dir_name_llsvp = sprintf('llsvp_%s', dir_name);
end
dir_out_base = '../../output/LLSVP';
dir_out_llsvp   = fullfile(dir_out_base, dir_name_llsvp);
mod_path_llsvp  = fullfile(dir_out_llsvp, mod_name);
            
% Create the output directory if it doesn't exist.
if ~exist(dir_out_llsvp, 'dir')
    
    mkdir(dir_out_llsvp);
    
end

%
if ~isnan(exagg)
    dir_name_ref = sprintf('prem_exagg_%s', dir_name);
else
    dir_name_ref = sprintf('prem_%s', dir_name);
end
dir_out_ref   = fullfile(dir_out_base, dir_name_ref);
mod_path_ref  = fullfile(dir_out_ref, mod_name);
            
% Create the output directory if it doesn't exist.
if ~exist(dir_out_ref, 'dir')
    
    mkdir(dir_out_ref);
    
end
            
% Save triangulated sphere as a .poly file.
trisurf2poly(mod_path_llsvp, pin, tin);
polyfile_llsvp = [mod_path_llsvp, '.poly'];
polyfile_ref = [mod_path_ref, '.poly'];

copyfile(polyfile_llsvp, polyfile_ref);

toc;

%% Create 3D mesh.

%l_max = (2.0^(1/2.0))*(3.0^(1/3.0))*(tet_max_vol^(1.0/3.0));
l_min = l_max/edge_length_ratio_mid_mantle_to_surf;

% Create the Tetgen mesh sizing function.
mesh_sizing_node_file   = [mod_path_llsvp, '.b.node'];
if isfile(mesh_sizing_node_file)

    fprintf('Tetgen mesh sizing node file already exists: %s\n', mesh_sizing_node_file)
    fprintf('Skipping creation of mesh sizing files.\n')
    
else

    %mesh_sizing(R1, R2, R3, l_min, l_max, mod_path_llsvp)
    %mesh_sizing(R1_adjusted, R2_adjusted, R3_adjusted, l_min, l_max, mod_path_llsvp)
    mesh_sizing_spheroid(R1, R2, R3, l_min, l_max, ellipticity_profile, mod_path_llsvp, dir_unit_spheroid, distmesh);
    
end

tic;

% Generate the mesh by calling tetgen from command line.
% This creates .ele, .neigh and .node files.
%unix([tetgen,' -pq1.5nYVFAa', num2str(tet_max_vol,'%f'),' ', mod_path,'.poly']);
% p     Take a piecewise linear complex (.poly file) and tetrahedralise it.
% q1.5/20.0 Improve quality of mesh so that each element has no radius/edge
%           ratio greater than 1.5, and dihedral angles no smaller than
%           20.0
% n         Output tetrahedra neighbours to .neigh file.
% Y         Do not modify input surface meshes.
% V         Verbose.
% F         Do not output .face or .edge files.
% A         Assign attributes to tetrahedra in different regions.
% a         Apply a maximum tetrahedron volume constraint.
%
unix([tetgen,' -A -pmq1.5/10.0nYCVFO5/7 -a ', num2str(tet_max_vol,'%f'),' ', mod_path_llsvp,'.poly']);

time = toc;

% Load tetgen output files (.ele, .node and .neigh).
% (See tetgen manual for more detail.)
% pout  Points of mesh.
% tout  Indices of points in each tetrahedron.
% at    'Attribute', value of 1, 2, or 3 indicating mantle, outer core
%       or inner core.
% neigh Neighbours. For each tetrahedron, gives the four indices of the
%       tetrahedra whose faces are opposite each of the four vertices.
%       Value of -1 indicates that face is on exterior.
% Can visualise mesh using tetramesh(tout, pout, 'FaceAlpha', 0.0).
[pout,tout,~,at,neigh] = read_mesh3d([mod_path_llsvp,'.1']);

% Write some additional output files.
fhed = [mod_path_llsvp, '.1_mesh.header'];
fele = [mod_path_llsvp, '.1_ele.dat'];
fngh = [mod_path_llsvp, '.1_neigh.dat'];
fnde = [mod_path_llsvp, '.1_node.dat'];
ftime = [mod_path_llsvp, '.1_time.dat'];

% The mesh.header file has the number of tetrahedra and number of nodes.
dh  = [size(tout,1) size(pout,1)];
fid = fopen(fhed,'w');
fprintf(fid,'%d %d',dh);

% The ele.dat file has the tetrahedra indices (in binary format).
fid = fopen(fele,'w');
fwrite(fid,tout','int');

% The neigh.dat file has the tetrahedra neighbours (in binary format).
fid = fopen(fngh,'w');
fwrite(fid,neigh','int');

% The node.dat file has the tetrahedra nodes (in binary format).
fid =fopen(fnde,'w');
fwrite(fid,pout','float64');

% The time file stores the time required to create tetrahedral mesh.
fid = fopen(ftime, 'w');
fprintf(fid, '%.1f', time);

fele_txt = [mod_path_llsvp, '.1.ele'];
fngh_txt = [mod_path_llsvp, '.1.neigh'];
fnde_txt = [mod_path_llsvp, '.1.node'];

fhed_ref = [mod_path_ref, '.1_mesh.header'];
fele_ref = [mod_path_ref, '.1_ele.dat'];
fngh_ref = [mod_path_ref, '.1_neigh.dat'];
fnde_ref = [mod_path_ref, '.1_node.dat'];
ftime_ref = [mod_path_ref, '.1_time.dat'];

fele_ref_txt = [mod_path_ref, '.1.ele'];
fngh_ref_txt = [mod_path_ref, '.1.neigh'];
fnde_ref_txt = [mod_path_ref, '.1.node'];

copyfile(fhed, fhed_ref)
copyfile(fele, fele_ref)
copyfile(fngh, fngh_ref)
copyfile(fnde, fnde_ref)
copyfile(ftime, ftime_ref)
copyfile(fele_txt, fele_ref_txt)
copyfile(fngh_txt, fngh_ref_txt)
copyfile(fnde_txt, fnde_ref_txt)

% Now, assign elastic parameters to the nodes.
make_vtk_file = 1; % Set this to 0 to skip creation of VTK file.
%mod_path_ref = 0; % Set this to 0 to skip creation of reference model.
%LLSVP_model(mod_path_llsvp, mod_path_ref, R2, pOrder, nlayer, at, MI, RD, make_vtk_file);
LLSVP_model_spheroid(mod_path_llsvp, mod_path_ref, pOrder, nlayer, at, MI, RD, ellipticity_profile, make_vtk_file);