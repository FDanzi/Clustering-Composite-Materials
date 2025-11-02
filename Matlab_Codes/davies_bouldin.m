% ========================================================================
%  Title: Cracking Composite Design: Trace-Driven Scaling for 
%         Cluster-Based Universal Design
%  Author: Francesco Danzi
%  Reference:
%     F. Danzi, "Cracking Composite Design: Trace-Driven Scaling for 
%     Cluster-Based Universal Design," Composite Structures, 2025.
%     DOI: [insert DOI once published]
%
% ------------------------------------------------------------------------
%  Description:
%     This script computes the Davies–Bouldin Index (DBI) to evaluate
%     the cohesiveness and separation of the identified clusters.
%
%  License:
%     MIT License (recommended for research sharing)
%
%  GitHub Repository:
%     https://github.com/FDanzi/Clustering-Composite-Materials
% ========================================================================



function dbi = davies_bouldin(X, labels)
%DAVIES_BOULDIN  Compute the Davies–Bouldin index
%
%   dbi = DAVIES_BOULDIN(X, labels)
%
%   Inputs:
%     X      – N×D data matrix (rows are points, columns are features)
%     labels – N×1 vector of integer cluster labels (1…k)
%
%   Output:
%     dbi    – scalar Davies–Bouldin index (lower is better)

    % Identify clusters
    C = unique(labels);
    k = numel(C);
    [N, ~] = size(X);
    
    % Preallocate
    centroids = zeros(k, size(X,2));
    S = zeros(k,1);       % within-cluster scatter (average distance to centroid)
    
    % Compute centroids and S_i
    for i = 1:k
        pts = X(labels==C(i), :);
        centroids(i,:) = mean(pts, 1);
        % average Euclidean distance of points in cluster i to centroid
        dists = sqrt(sum((pts - centroids(i,:)).^2, 2));
        S(i) = mean(dists);
    end
    
    % Compute pairwise centroid distances M_ij
    M = squareform( pdist(centroids, 'euclidean') );
    
    % Compute R_ij = (S_i + S_j) / M_ij for i ≠ j
    R = zeros(k);
    for i = 1:k
      for j = 1:k
        if i ~= j
          R(i,j) = (S(i) + S(j)) / M(i,j);
        end
      end
    end
    
    % For each i, find the maximum R_ij over j ≠ i
    R_max = max(R, [], 2);
    
    % Davies–Bouldin index is the average of those maxima
    dbi = mean(R_max);
end
