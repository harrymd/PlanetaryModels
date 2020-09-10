function [volume] = volume_of_triangulation(indices, vertices)
    %%
    % From https://www.mathworks.com/matlabcentral/fileexchange/26982-volume-of-a-surface-triangulation
    %
    % vertices     (nPoints x 3)       verticesertices of triangulation.
    % indices     (nTriangles x 3)    Indices of triangulation.
    
    % 'z' coordinate of triangle centers.
    FaceCentroidZ = ( vertices(indices(:, 1), 3) + vertices(indices(:, 2), 3) + vertices(indices(:, 3), 3) ) /3;
    % Face normal vectors, with length equal to triangle area.
    FNdA = cross( (vertices(indices(:, 2), :) - vertices(indices(:, 1), :)), ...
    (vertices(indices(:, 3), :) - vertices(indices(:, 2), :)) , 2 ) / 2;
    % verticesolume from divergence theorem (using vector field along z).
    volume = FaceCentroidZ' * FNdA(:, 3);
end

