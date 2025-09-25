    %% Run to get the filename.gml file into the datasets folder as filename.mat

    filename='football';
    S = fileread(strcat('./datasets\',filename,'.gml'));
    nodes = regexp(S, 'node.*?id (?<id>\d+).*?label\s*"(?<label>[^"]*)"', 'names');
    edges = regexp(S, 'edge.*?source\s*(?<source>\d+).*?target\s*(?<target>\d+)', 'names');
    all_ids = {nodes.id};
    all_names = {nodes.label};
    all_sources = {edges.source};
    all_targets = {edges.target};
    [source_found, s] = ismember(all_sources, all_ids);
    nfidx = find(~source_found);
    if ~isempty(nfidx)
       error('Source ids not found in node list, starting from "%s"', edges(nfidx(1).source));
    end
    [target_found, t] = ismember(all_targets, all_ids);
    nfidx = find(~target_found);
    if ~isempty(nfidx)
       error('Target ids not found in node list, starting from "%s"', edges(nfidx(1).target));
    end
    EdgeTable = table([s.', t.'], ones(length(s),1), 'VariableNames', {'EndNodes' 'Weight'});
    NodeTable = table(all_names.', 'VariableNames',{'Name'});
    G = graph(EdgeTable,NodeTable);
    A=adjacency(G);
    save(strcat('./datasets\',filename,'.mat'),'A');