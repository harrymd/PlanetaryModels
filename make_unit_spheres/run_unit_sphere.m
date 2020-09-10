function [] = run_unit_sphere(vol_max_tet, out_dir, distmesh_path, radius)
%matlab -nojvm -r "tet_max_vol=8.0E10; radius=6371.0; out_dir = './output';  distmesh_path = '../PlanetaryModels/distmesh'; try; run_unit_sphere(tet_max_vol, out_dir, distmesh_path, radius); catch; exit; end; exit"
    [initial_spacing]  = tet_max_vol_to_distmesh_spacing(vol_max_tet, radius);
    [file_name]         = unit_sphere(initial_spacing, out_dir, distmesh_path);
end