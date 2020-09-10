function [file_name] = unit_spheroid(initial_spacing, epsilon, out_dir, distmesh_path)
    %% Build an approximately uniform mesh on a "unit spheroid".
    % A "unit spheroid" is a unixial ellipsoid with two principal axes (the
    % x and y axes) of length 1.0, and a third (shorter) axis (the z axis)
    % of length (1.0 - epsilon), where epsilon is the ellipticity.
    %
    % Input
    % initial_spacing   This is an input to the distmesh function which
    %                   controls how fine the grid will be. Use the
    %                   function tet_max_vol_to_distmesh_spacing() to
    %                   calculate the appropriate value for a given maximum
    %                   tetrahedron volume.
    % epsilon           The ellipticity.
    % out_dir           The output directory.
    % distmesh_path     Path to the distmesh function.
    
    % Check if the mesh already exists.
    name_str = sprintf('%12.6E_%10.8f', initial_spacing, epsilon);
    file_name = sprintf('unit_spheroid_%s.mat', name_str);
    file_name = fullfile(out_dir, file_name);
    
    prefix = sprintf('%sunit_spheroid_%s', out_dir, name_str);
    file_name_nodes = sprintf('%s_nodes.txt', prefix);
    file_name_elements = sprintf('%s_elements.txt', prefix);
    
    if isfile(file_name)
        
        fprintf('Unit spheroid mesh already exists: %s\n', file_name)
        fprintf('Skipping mesh creation.\n')
        
    else
        
        fprintf('Building unit spheroid mesh with spacing parameter %12.6E and ellipticity %10.8f\n', initial_spacing, epsilon)

        % Add the distmesh package.
        addpath(distmesh_path)

        % Load the distance function for a unit spheroid at the origin.
        % This is supposed to give the signed distance from a point to the surface of the spheroid, such that result is negative inside the spheroid, 0 on the boundary and positive
        % outside the spheroid. However, this example (based on the
        % official distmesh documentation), seems to be always positive and
        % not equal to the distance to the spheroid.
        fd=@(p) (p(:,1).^2.0) + (p(:,2).^2.0) + (p(:,3)/(1.0 - epsilon)).^2 - 1.0;

        % Start the timer.
        tic

        % Use distmesh to calculate the mesh on the surface of the spheroid.
        % p     points
        % t     triangulation
        [p,t] = distmeshsurface(                    ...
                    fd,                             ...
                    @huniform,                      ...
                    initial_spacing,                ...
                    1.1*[-1, -1, -1; 1, 1, 1]);
        
        vol = abs(volume_of_triangulation(t, p));
                
        time = toc;
        % Save the variables to a Matlab file.
        save(file_name, 'fd', 'p', 'initial_spacing', 't', 'time', 'vol');
        
    end
    
end