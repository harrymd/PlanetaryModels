function [file_name] = unit_sphere(initial_spacing, out_dir, distmesh_path)
    %% Build an approximately uniform mesh on a unit sphere.
    %
    % Input
    % initial_spacing   This is an input to the distmesh function which
    %                   controls how fine the grid will be. Use the
    %                   function tet_max_vol_to_distmesh_spacing() to
    %                   calculate the appropriate value for a given maximum
    %                   tetrahedron volume.
    % out_dir           The output directory.
    % distmesh_path     Path to the distmesh function.
    
    % Check if the mesh already exists.
    file_name = sprintf('unit_sphere_%12.6E.mat', initial_spacing);
    file_name = fullfile(out_dir, file_name);
    
    prefix = sprintf('%sunit_sphere_%12.6E', out_dir, initial_spacing);
    file_name_nodes = sprintf('%s_nodes.txt', prefix);
    file_name_elements = sprintf('%s_elements.txt', prefix);
    
    if isfile(file_name)
        
        fprintf('Unit sphere mesh already exists: %s\n', file_name)
        fprintf('Skipping mesh creation.\n')
        
    else
        
        fprintf('Building unit sphere mesh with spacing parameter %12.6E\n', initial_spacing)

        % Add the distmesh package.
        addpath(distmesh_path)

        % Load the distance function for a unit sphere at the origin.
        % This gives the signed distance from a point to the surface of the sphere.
        % The result is negative inside the sphere, 0 on the boundary and positive
        % outside the sphere.
        fd = @(p) dsphere(p, 0, 0, 0, 1);

        % Start the timer.
        tic

        % Use distmesh to calculate the mesh on the surface of the sphere.
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
        
        %file_ID_nodes = fopen(file_name_nodes, 'w');
        %fprintf(file_ID_nodes, '%>19.12e %>19.12e %>19.12e\n', p');
        %fclose(file_ID_nodes);
        
        %file_ID_elements = fopen(file_name_elements, 'w');
        %fprintf(file_ID_elements, '%12d %12d %12d\n', t');
        %fclose(file_ID_elements);
        
        %triangulation_orient(prefix)

    end
    
end
% for reproducing

% surf_d = 0.3;   p 198;   t 392
% surf_d = 0.2;   p 480;   t 956
% surf_d = 0.1;   p 1806;  t 3608
% surf_d = 0.08;  p 3k;    t 6k
% surf_d = 0.05;  p 7k;    t 15k
% surf_d = 0.03;  p 21k;   t 42k
% surf_d = 0.02;  p 47k;   t 94k
% surf_d = 0.015; p 84k;   t 167k
% surf_d = 0.01;  p 188k;  t 377k
% surf_d = 0.008; p 294k;  t 589k
% surf_d = 0.006; p 523k;  t 1047k?
% surf_d = 0.004; p 1177k; t 2353k time:361294.354594
% surf_d = 0.0035; p 1538k; t 3077k time: 771746.288223


