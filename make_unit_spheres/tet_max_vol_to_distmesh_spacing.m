function [distmesh_spacing] = tet_max_vol_to_distmesh_spacing(vol_max_tet, radius)

    % Estimate the number of triangles on the free surface.
    % First option: Based on ratio of surface areas between the sphere and
    % a tetrahedral element.
    %(Seems to give too few tetrahedra,
    % so we multiply by 30 to give better results.)
    area_sphere         = 4.0*pi*(radius^2.0);
    area_max_tet_face   = 8.0*(3.0^(-7.0/6.0))*(vol_max_tet^(2.0/3.0));
    n_triangles         = 30.0*round(area_sphere/area_max_tet_face);
    % Second option: Empirical.
    % n_tri_surf = 2.0E8/sqrt(tet_max_vol);
    % n_tri_surf  = round(n_tri_surf);
    
    % Calculate the edge-length parameter (d) for the distmesh function. This
    % parameter controls the number of triangles in the surface meshing:
    % approximately, n = k/(d**2) where k = 37.67795182. This approximation is
    % accurate to within 10% for small meshes, and the accuracy improves for
    % larger ones.
    sqrtk               = 6.138236866;
    distmesh_spacing    = sqrtk/sqrt(n_triangles);
    % A bug in distmeshsurface can cause it to fail if the initial grid is not
    % an exact multiple of the d parameter (initial spacing). Therefore, 
    % we change the d parameter so it is equal to the nearest integer divisor.
    d_grid_width     = 2.2;
    distmesh_spacing = d_grid_width/(round(d_grid_width/distmesh_spacing));

end