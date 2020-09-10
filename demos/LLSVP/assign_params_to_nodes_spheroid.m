function [vp0, vs0, rho0] = assign_params_to_nodes(         ...
                                pNp, nlayer, x, y, z,       ...
                                at, MI, RD,                 ...
                                ellipticity_profile,        ...
                                r_max_llsvp,                ...
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
    
    % Create substitute for Matlab's built-in legendre() function, which
    % calculates more output than necessary.
    % legendreP2 Legendre polynomial of degree two.
    legendreP2 = @(x) ((3.0*(x.^2)) - 1.0)/2.0;

    % Assign values at nodes using interpolation of input model.
    for ati = 1:nlayer
        
        % Find indices of nodes in the given layer (inner core, outer core
        % or mantle).
        tid = find(at==ati); 

        % Find radial coordinate of nodes.
        rprime_vd = sqrt(x(:,tid).^2+y(:,tid).^2+z(:,tid).^2);
        
        % Find the polar angle of the nodes.
        theta = acos(z(:, tid)./rprime_vd);
        P2costheta = legendreP2(cos(theta));
        
        % Convert radial coordinate (rprime) to original (unflattened)
        % radial coordinate (r). This involves solving Dahlen and Tromp
        % (1998) eq. 14.4 for r. We make the approximation that 
        % eps(r) ~ eps(rprime)
        % which leads to a simple rearrangement.
        % First, find ellipticity (applying the approximation).
        % Note that some points fall outside the profile due to flattening.
        % These are given the outer surface ellipticity.
        eps_vd = interp1(ellipticity_profile(1, :), ellipticity_profile(2, :), rprime_vd, 'linear', ellipticity_profile(2, end));
        % Second, calculate the original radial coordinates.
        r_vd = rprime_vd./(1.0 - ((2.0/3.0)*(eps_vd.*P2costheta)));
        
        % Can calculate a second approximation of the ellipticity and
        % therefore the radius.
        eps_vd_updated = interp1(ellipticity_profile(1, :), ellipticity_profile(2, :), r_vd, 'linear', ellipticity_profile(2, end));
        r_vd_updated = rprime_vd./(1.0 - ((2.0/3.0)*(eps_vd_updated.*P2costheta)));
        disp(max(abs(r_vd(:) - r_vd_updated(:))));
        eps_vd = eps_vd_updated;
        r_vd = r_vd_updated;
        
        %
        %eps_vd_updated = interp1(ellipticity_profile(1, :), ellipticity_profile(2, :), r_vd, 'linear', ellipticity_profile(2, end));
        %r_vd_updated = rprime_vd./(1.0 - ((2.0/3.0)*(eps_vd_updated.*P2costheta)));
        %disp(max(abs(r_vd(:) - r_vd_updated(:))));
        
        % The labelling of attributes by TetGen is inconsistent, so we
        % have to work out which attribute corresponds to which layer.
        r_vd_flat = reshape(r_vd, numel(r_vd), 1);
        r_vd_min = min(r_vd_flat);
        r_vd_max = max(r_vd_flat);
        r_vd_mid = (r_vd_min + r_vd_max)/2.0;
        if r_vd_mid < RD(3, 1)
            
            i = 3;
           
        elseif r_vd_mid < RD(2, 1)
                
            i = 2;
            
        elseif r_vd_mid < RD(1, 1)
            
            i = 1;
         
        else
            
            i = NaN;
            
        end

        % Extract the portion of the radial model from this layer.
        tmpr   = MI(RD(i,2) + 1 : RD(i+1,2), 1);
        tmprho = MI(RD(i,2) + 1 : RD(i+1,2), 2);
        tmpvp  = MI(RD(i,2) + 1 : RD(i+1,2), 3);
        tmpvs  = MI(RD(i,2) + 1 : RD(i+1,2), 4);
        
        %fprintf('%1d %7.1f %7.1f %7.1f %7.1f%', i, min(tmpr), max(tmpr), min(min(r_vd(:, :))), max(max(r_vd(:, :))))
        %disp('\n')
        
        % Find Vp by interpolation.
        % The flag 'pchip' indicates shape-preserving cubic spline
        % interpolation.
        vptmp       = interp1(tmpr,tmpvp,r_vd,'pchip');

        % Find Vs by interpolation.
        vstmp       = interp1(tmpr,tmpvs,r_vd,'pchip');

        % Find rho by interpolation.
        rhotmp      = interp1(tmpr,tmprho,r_vd,'pchip');
        
        %disp((min(reshape(vptmp, numel(vptmp), 1))))
        %disp((max(reshape(vptmp, numel(vptmp), 1))))

        % If requested, add the LLSVP to the mantle.
        if ~isnan(lon_llsvp)
            
            % Add the LLSVP to the mantle.
            if i == 1

                % Find longitude and latitude of the nodes.
                [lon, lat , ~] = cart2sph(x(:,tid), y(:,tid), z(:,tid));
                
                % Convert from geocentric to geographic coordinates.
                % (Dahlen and Tromp, 1998, section 14.1.5.)
                colat = pi/2.0 - lat;
                tan_colat_prime = tan(colat)./(1.0 + 2.0*eps_vd);
                colat_prime = atan(tan_colat_prime);
                lat_prime = pi/2.0 - colat_prime;
                
                % Fix angle wrap-around convention.
                lat_prime_out_of_bounds = (lat_prime > pi/2.0);
                lat_prime(lat_prime_out_of_bounds) = lat_prime(lat_prime_out_of_bounds) - pi;
                
                % Convert to degrees.
                lon = rad2deg(lon);
                lat_prime = rad2deg(lat_prime);

                % Check if the nodes lie within the outline of the LLSVP.
                in_llsvp_poly   = inpolygon(lon, lat_prime, lon_llsvp, lat_llsvp);

                % Check if the nodes lie within the radial extent of the LLSVP.
                % (Note: Use undeformed radial coordinate).
                in_llsvp_rad    = ((r_vd >= RD(2, 1)) & (r_vd <= r_max_llsvp));

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