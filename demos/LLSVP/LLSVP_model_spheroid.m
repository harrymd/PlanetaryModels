function [] = LLSVP_model_spheroid(mod_path, mod_path_ref, pOrder, nlayer, at, MI, RD, ellipticity_profile, make_vtk)

    tic;
    
    % Load the outline of the LLSVP.
    file_llsvp  = 'llsvp_african_outline.txt';
    data_llsvp  = load(file_llsvp);
    lon_llsvp   = data_llsvp(:, 1);
    lat_llsvp   = data_llsvp(:, 2);

    % Define the anomalies in Vp, Vs and rho (as a fraction of reference value.
    dVp_over_Vp =   -0.016;
    dVs_over_Vs =   -0.040;
    drho_over_rho = +0.010;

    % Define the height of the LLSVP (km).
    d_r_llsvp   = 400.0;
    r_max_llsvp = RD(2, 1) + d_r_llsvp;

    % Declare file names.
    fname   = [mod_path,'.1']; 
    %
    fmid    = ['_pod_',int2str(pOrder),'_'];
    ftail   = 'true.dat';
    %
    fvp     = [fname,'_vp', fmid,ftail];
    fvs     = [fname,'_vs', fmid,ftail];
    frho    = [fname,'_rho',fmid,ftail];
    fvtk    = [fname,fmid,'model.vtk'];

    % pNp Number of nodes per tetrahedra (4 if pOrder = 1,
    % 10 if pOrder = 2).
    pNp = (pOrder+1)*(pOrder+2)*(pOrder+3)/6;
    % (nele)            Number of elements.
    % tet [nele*8, 4]   Indices for sub-elements (?).
    % x, y, z [pNp, nele] CCoordinates of each node.
    [x,y,z,tet] = construct(fname,pOrder);

    [vp0, vs0, rho0]                = assign_params_to_nodes_spheroid(      ...
                                        pNp, nlayer,                        ...
                                        x, y, z,                            ...
                                        at, MI, RD,                         ...
                                        ellipticity_profile,                ...
                                        r_max_llsvp,                        ...
                                        lon_llsvp, lat_llsvp,               ...
                                        dVp_over_Vp, dVs_over_Vs,           ...
                                        drho_over_rho);

    % Build a model with no LLSVP.                    
    [vp0_ref, vs0_ref, rho0_ref]    = assign_params_to_nodes_spheroid(      ...
                                        pNp, nlayer,                        ...
                                        x, y, z,                            ...
                                        at, MI, RD,                         ...
                                        ellipticity_profile,                ...
                                        NaN,                                ...
                                        NaN, NaN,                       	...
                                        NaN, NaN,                        	...
                                        NaN);

    % Build arrays of anomalies (for convenience).
    dvp0    = vp0 - vp0_ref;
    dvs0    = vs0 - vs0_ref;
    drho0   = rho0 - rho0_ref;

    % Build arrays of fractional anomalies (for convenience).
    dvp0_over_vp0   = dvp0./vp0_ref;
    dvs0_over_vs0   = dvs0./vs0_ref;
    drho0_over_rho0 = drho0./rho0_ref;

    % Correct for division-by-zero error.
    dvs0_over_vs0(isnan(dvs0_over_vs0)) = 0.0;

    % Write model to files.
    accry   = 'float64';
    %
    fid     = fopen(fvp,'w');
    fwrite(fid, vp0(:), accry);
    fclose(fid);
    %
    fid     = fopen(fvs,'w');
    fwrite(fid, vs0(:), accry);
    fclose(fid);
    %
    fid     = fopen(frho,'w');
    fwrite(fid, rho0(:), accry);
    fclose(fid);
    
    time = toc;
    time_file = [mod_path, '_param_assign_time.dat'];
    fid = fopen(time_file, 'w');
    fprintf(fid, '%.1f', time);
    fclose(fid);
    
    if mod_path_ref ~= 0
        
         fname   = [mod_path_ref,'.1']; 
 
         fvp     = [fname,'_vp', fmid,ftail];
         fvs     = [fname,'_vs', fmid,ftail];
         frho    = [fname,'_rho',fmid,ftail];
         
         fid     = fopen(fvp,'w');
         fwrite(fid, vp0_ref(:), accry);
         fclose(fid);
         %
         fid     = fopen(fvs,'w');
         fwrite(fid, vs0_ref(:), accry);
         fclose(fid);
         %
         fid     = fopen(frho,'w');
         fwrite(fid, rho0_ref(:), accry);
         fclose(fid);
         
     end

    %% Write VTK file.
    if make_vtk == 1
        
        % Re-order (I'm not sure why...)
        vp  = vp0(tet');
        vs  = vs0(tet');
        rho = rho0(tet');

        % Names.
        filename = fvtk;
        data_title = 'model';

        % Create data structure.
        % .type field must be 'scalar' or 'vector'.
        % If vector, the number of components must equal the dimension of the mesh
        data_struct(1).type = 'scalar';
        data_struct(1).name = 'Vp';
        data_struct(1).data = vp(:);

        data_struct(2).type = 'scalar';
        data_struct(2).name = 'Vs';
        data_struct(2).data = vs(:);

        data_struct(3).type = 'scalar';
        data_struct(3).name = 'Density';
        data_struct(3).data = rho(:);

        % Additional fields (for plotting).
        dvp             = dvp0(tet');
        dvs             = dvs0(tet');
        drho            = drho0(tet');
        %
        dvp_over_vp     = dvp0_over_vp0(tet');
        dvs_over_vs     = dvs0_over_vs0(tet');
        drho_over_rho   = drho0_over_rho0(tet');
        %
        data_struct(4).type = 'scalar';
        data_struct(4).name = 'dVs';
        data_struct(4).data = dvs(:);
        %
        data_struct(5).type = 'scalar';
        data_struct(5).name = 'dVp';
        data_struct(5).data = dvp(:);
        %
        data_struct(6).type = 'scalar';
        data_struct(6).name = 'drho';
        data_struct(6).data = drho(:);
        %
        data_struct(7).type = 'scalar';
        data_struct(7).name = 'dVp_Vp';
        data_struct(7).data = dvp_over_vp(:);
        %
        data_struct(8).type = 'scalar';
        data_struct(8).name = 'dVs_Vs';
        data_struct(8).data = dvs_over_vs(:);
        %
        data_struct(9).type = 'scalar';
        data_struct(9).name = 'drho_rho';
        data_struct(9).data = drho_over_rho(:);

        % Set to true to flip data.
        flipped = false;

        % Get re-ordered list of points.
        tnew        = reshape(1:size(tet,1)*size(tet,2),size(tet,2),size(tet,1));
        psiz        = max(tet(:));
        pnew0       = zeros(psiz,3);
        pnew0(:,1)  = x(:); 
        pnew0(:,2)  = y(:); 
        pnew0(:,3)  = z(:);  
        pnew        = pnew0(tet',:);

        % Write the VTK file.
        stat = vtk_write_tetrahedral_grid_and_data(filename,data_title,pnew,...
            tnew',data_struct,flipped);

    end

end