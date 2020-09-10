function [] = mesh_sizing_spheroid(r1, r2, r3, l_min, l_max, ellipticity_profile, file_root, dir_unit_spheroid, distmesh)
%% Builds a mesh sizing function for a spheroid.
% eps Ellipticity.
% r1 > r2 > r3
    
    tic;

    %% Create the nodes and the triangulation.
    
    r_gap_max   = max([abs(r1 - r2), abs(r2 - r3)]);
    l_range     = l_max - l_min;
    dl_dr       = l_range/(r_gap_max/2.0);
    
    r_list = [0.0, r3/2.0, r3, (r2 + r3)/2.0, r2, (r2 + r1)/2.0, r1, r1*1.05];
    l_list = [0.0, dl_dr*((r3 - 0.0)/2.0), 0.0, dl_dr*(r2 - r3)/2.0, 0.0, dl_dr*(r1 - r2)/2.0, 0.0, 0.0];
    l_list = l_min + l_list;
    
    eps_list = interp1(ellipticity_profile(1, :), ellipticity_profile(2, :), r_list, 'linear', ellipticity_profile(2, end));
    
    % Set the spheroid grid side.
    d_surf = 0.22; % Gives about 400 nodes on the spheroid.
    
    % Loop over array.
    n = length(r_list);
    % Initialise output lists.
    % The first item is a single point at the centre of the planet.
    p = [0.0, 0.0, 0.0];
    l = l_min;
    % Note loop starts from 2 (centre of planet was i = 1).
    for i = 2 : n
        
        % Generate the unit spheroid.
        unit_spheroid_file_i   = unit_spheroid(d_surf, eps_list(i), dir_unit_spheroid, distmesh);
        pi_struct = load(unit_spheroid_file_i, 'p');
        pi = pi_struct.p;
        
        % Scale.
        r_equator_i = r_list(i)*(1.0 + (eps_list(i)/3.0));
        pi         = r_equator_i*pi; 

        % Set a constant edge length for each point on this spheroid.
        [ni, ~] = size(pi);
        li = l_list(i)*ones(ni, 1)';
        
        % Store results from this iteration.
        p = cat(1, p, pi);
        l = [l, li];
        
    end
    
    clear pi li
    
    %p = [x; y; z]';
    %[p, ia, ~] = unique(p, 'first', 'rows');
    %l = l(ia);

    tri = delaunayTriangulation(p);
    
    % Number of nodes and tetrahedra.
    n = length(l);
    nt = tri.size(1);
    
%     % Check interpolation.
%     nr = 100000;
%     pr = 2.0*r1*(rand(3, nr) - 0.5)';
%     prr = sqrt(pr(:, 1).^2.0 + pr(:, 2).^2.0 + pr(:, 3).^2.0);
%     pr = pr(prr < r1, :);
%     prr = prr(prr < r1);
%     
%     [ti,bc] = pointLocation(tri,pr);
%     triVals = l(tri(ti,:));
%     
%     Vq = dot(bc',triVals')';
%     
%     scatter(prr, Vq)
    

    %% Save output to files.
    
    node_file   = [file_root, '.b.node'];
    ele_file    = [file_root, '.b.ele'];
    mtr_file    = [file_root, '.b.mtr'];
    time_file   = [file_root, '.b.time'];
    
    % Node file.
    fid = fopen(node_file, 'w');
    
    fprintf(fid, '%12i %12i %12i %12i\n', n, 3, 0, 0);
    
    for i = 1 : n
       
        fprintf(fid, '%12i %12.5e %12.5e %12.5e\n', i, tri.Points(i, 1), tri.Points(i, 2), tri.Points(i, 3));
        
    end
    
    fclose(fid);
    
    % Element file.
    
    fid = fopen(ele_file, 'w');
    
    fprintf(fid, '%12i %12i %12i\n', nt, 4, 0);
    
    for i = 1 : nt
       
        fprintf(fid, '%12i %12i %12i %12i %12i\n', i, tri(i, 1), tri(i, 2), tri(i, 3), tri(i, 4));
        
    end
    
    fclose(fid);
    
    % Mesh file.
    fid = fopen(mtr_file, 'w');
    
    fprintf(fid, '%12i %12i\n', n, 1);
    
    for i = 1 : n
       
        fprintf(fid, '%12.5e\n', l(i));
        
    end
    
    fclose(fid);
    
    time = toc;
    fid = fopen(time_file, 'w');
    
    fprintf(fid, '%.1f', time);

end