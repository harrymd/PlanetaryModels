function [] = mesh_sizing(r1, r2, r3, l_min, l_max, file_root)
%% Builds a mesh sizing function.
% r1 > r2 > r3
    
    tic;

    %% Create the nodes and the triangulation.
    
    r_gap_max   = max([abs(r1 - r2), abs(r2 - r3)]);
    l_range     = l_max - l_min;
    dl_dr       = l_range/(r_gap_max/2.0);
    
    r_list = [0.0, r3/2.0, r3, (r2 + r3)/2.0, r2, (r2 + r1)/2.0, r1, r1*1.05];
    l_list = [0.0, dl_dr*((r3 - 0.0)/2.0), 0.0, dl_dr*(r2 - r3)/2.0, 0.0, dl_dr*(r1 - r2)/2.0, 0.0, 0.0];
    l_list = l_min + l_list;
    
    n = length(r_list);
    x = 0.0;
    y = 0.0;
    z = 0.0;
    l = l_min;
    n_p = 20;
    for i = 2 : n
        
        % % This is too slow.
        % n_p = ceil(2.0*pi*r_list(i)/l_min);
        
        [xi, yi, zi] = sphere(n_p);
        
        % Flatten and scale.
        xi = r_list(i)*xi(:)';
        yi = r_list(i)*yi(:)';
        zi = r_list(i)*zi(:)';
        
        ni = length(xi);
        li = l_list(i)*ones(ni, 1)';
        
        x = [x, xi];
        y = [y, yi];
        z = [z, zi];
        l = [l, li];
        
    end
    
    clear xi yi zi li
    
    p = [x; y; z]';
    [p, ia, ~] = unique(p, 'first', 'rows');
    l = l(ia);

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