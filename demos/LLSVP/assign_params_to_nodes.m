function [vp0, vs0, rho0] = assign_params_to_nodes(         ...
                                pNp, nlayer, x, y, z,       ...
                                at, MI, RD, r_max_llsvp,    ...
                                lon_llsvp, lat_llsvp,       ...
                                dVp_over_Vp, dVs_over_Vs,   ...
                                drho_over_rho)
    %% Builds a planetary model with an LLSVP.
    % First, builds a spherically-symmetric model by interpolating a radial
    % reference model at the nodes.
    % Then, if requested, a perturbation is added to all nodes lying within
    % the specified region.

    vp0  = zeros(pNp,size(x,2)); 
    vs0  = zeros(pNp,size(x,2)); 
    rho0 = zeros(pNp,size(x,2)); 

    % Assign values at nodes using interpolation of input model.
    i = 0;
    %for i = 1:nlayer
    %for j = [2, 1, 3]
    for ati = 1:nlayer
        
        % Find indices of nodes in the given layer (inner core, outer core
        % or mantle).
        tid = find(at==ati); 

        % Find radial coordinate of nodes.
        rvd = sqrt(x(:,tid).^2+y(:,tid).^2+z(:,tid).^2);
        
        % The labelling of attributes by TetGen is inconsistent, so we
        % have to work out which attribute corresponds to which layer.
        rvd_flat = reshape(rvd, numel(rvd), 1);
        rvd_min = min(rvd_flat);
        rvd_max = max(rvd_flat);
        rvd_mid = (rvd_min + rvd_max)/2.0;
        if rvd_mid < RD(3, 1)
            
            i = 3;
           
        elseif rvd_mid < RD(2, 1)
                
            i = 2;
            
        elseif rvd_mid < RD(1, 1)
            
            i = 1;
         
        else
            
            i = NaN;
            
        end

        % Extract the portion of the radial model from this layer.
        tmpr   = MI(RD(i,2) + 1 : RD(i+1,2), 1);
        tmprho = MI(RD(i,2) + 1 : RD(i+1,2), 2);
        tmpvp  = MI(RD(i,2) + 1 : RD(i+1,2), 3);
        tmpvs  = MI(RD(i,2) + 1 : RD(i+1,2), 4);
        
        %fprintf('%1d %7.1f %7.1f %7.1f %7.1f%', i, min(tmpr), max(tmpr), min(min(rvd(:, :))), max(max(rvd(:, :))))
        %disp('\n')
        
        % Find Vp by interpolation.
        % The flag 'pchip' indicates shape-preserving cubic spline
        % interpolation.
        vptmp       = interp1(tmpr,tmpvp,rvd,'pchip');

        % Find Vs by interpolation.
        vstmp       = interp1(tmpr,tmpvs,rvd,'pchip');

        % Find rho by interpolation.
        rhotmp      = interp1(tmpr,tmprho,rvd,'pchip');
        
        %disp((min(reshape(vptmp, numel(vptmp), 1))))
        %disp((max(reshape(vptmp, numel(vptmp), 1))))

        
        % If requested, add the LLSVP to the mantle.
        if ~isnan(lon_llsvp)
            
            % Add the LLSVP to the mantle.
            if i == 1

                % Find longitude and latitude of the nodes.
                [lon, lat , ~] = cart2sph(x(:,tid), y(:,tid), z(:,tid));
                lon = rad2deg(lon);
                lat = rad2deg(lat);

                % Check if the nodes lie within the outline of the LLSVP.
                in_llsvp_poly   = inpolygon(lon, lat, lon_llsvp, lat_llsvp);

                % Check if the nodes lie within the radial extent of the LLSVP.
                in_llsvp_rad    = ((rvd >= RD(2, 1)) & (rvd <= r_max_llsvp));

                % Find the nodes which satisfy both conditions.
                in_llsvp        = (in_llsvp_poly & in_llsvp_rad);

                % Apply the perturbation to the parameters of the nodes inside
                % the LLSVP.
                vptmp(in_llsvp)     = vptmp(in_llsvp)*(1.0 + dVp_over_Vp);
                vstmp(in_llsvp)     = vstmp(in_llsvp)*(1.0 + dVs_over_Vs);
                rhotmp(in_llsvp)    = rhotmp(in_llsvp)*(1.0 + drho_over_rho);

            end
        
        end

        % Put into final output arrays.
        vp0(:,tid)  = reshape(vptmp,pNp,length(tid));
        vs0(:,tid)  = reshape(vstmp,pNp,length(tid));
        rho0(:,tid) = reshape(rhotmp,pNp,length(tid));

    end

end