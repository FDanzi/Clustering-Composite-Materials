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
%     This script computes the Silhouette Index to assess the 
%     cohesiveness and separation of the identified clusters.
%
%  License:
%     MIT License (recommended for research sharing)
%
%  GitHub Repository:
%     https://github.com/FDanzi/Clustering-Composite-Materials
% ========================================================================



function [s, s_avg, s_cluster] = my_silhouette(X, labels)
%MY_SILHOUETTE  Compute silhouette scores for clustered data
%
%   [s, s_avg, s_cluster] = MY_SILHOUETTE(X, labels)
%
%   Inputs:
%     X         – N×D data matrix (rows are points, columns are features)
%     labels    – N×1 vector of integer cluster labels (1…k)
%
%   Outputs:
%     s         – N×1 vector of silhouette scores s(i) ∈ [–1,1]
%     s_avg     – scalar average silhouette over all points
%     s_cluster – k×1 vector of average silhouette per cluster

    % Number of points
    N = size(X,1);
    % Unique clusters
    C = unique(labels);
    k = numel(C);
    
    % Precompute full pairwise distance matrix
    D = squareform( pdist(X, 'euclidean') );
    
    % Initialize per‐point silhouette
    s = zeros(N,1);
    
    for i = 1:N
        ci = labels(i);
        % Indices of points in the same cluster (excluding i)
        same = find(labels == ci);
        same(same == i) = [];
        
        % 1) a(i): average distance to all other points in its cluster
        if isempty(same)
            a_i = 0;
        else
            a_i = mean( D(i, same) );
        end
        
        % 2) b(i): minimum average distance to any other cluster
        b_i = inf;
        for cj = setdiff(C, ci)'
            group_j = find(labels == cj);
            d_ij = mean( D(i, group_j) );
            if d_ij < b_i
                b_i = d_ij;
            end
        end
        
        % 3) silhouette score
        if isempty(same)
        s(i) = 0;    % force singleton silhouettes to zero
        else
        s(i) = (b_i - a_i)/max(a_i, b_i);
        end
    end
    
    % Overall average silhouette
    s_avg = mean(s);
    
    % Average silhouette per cluster
    s_cluster = zeros(k,1);
    for j = 1:k
        members = (labels == C(j));
        s_cluster(j) = mean( s(members) );
    end
end