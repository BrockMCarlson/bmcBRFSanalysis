% Example correlation matrix (replace with your data)
correlation_matrix = rand(32, 32);

% Apply a correlation threshold
threshold = 0.7;
correlation_matrix(correlation_matrix < threshold) = 0;

% Create a directed graph from the correlation matrix
G = graph(correlation_matrix, 'upper', 'OmitSelfLoops');

% Calculate network metrics
degree = centrality(G, 'degree');
clustering_coefficient = clusteringcoef(G);
betweenness_centrality = centrality(G, 'betweenness');

% Rest of the code remains the same...


% Perform community detection
communities = community_louvain(G);

% Visualize the network
figure;
h = plot(G, 'Layout', 'force', 'NodeCData', communities, 'NodeColor', 'flat');
colormap(jet);
colorbar;
title('Network Visualization');

% Display network metrics
disp('Degree:');
disp(degree);
disp('Clustering Coefficient:');
disp(clustering_coefficient);
disp('Betweenness Centrality:');
disp(betweenness_centrality);
disp('Community Assignments:');
disp(communities);
