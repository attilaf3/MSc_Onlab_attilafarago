function [lodf] = lodf(n2sPTDF, slack, G, Edge_table)

    % slack oszlop visszahelyezése a ptdf-be
    slack_id = find(G.Nodes.Name == slack);
    n2sPTDF = [n2sPTDF(:,1:slack_id-1), zeros(height(G.Edges),1), n2sPTDF(:,slack_id:end)];

    lodf = zeros(height(G.Edges));
    [~, loc] = ismember(Edge_table.Name, G.Edges.Name);
    end_nodes = string(G.Edges.EndNodes(loc,:));
    [~,node_s_id] = ismember(end_nodes(:,1), G.Nodes.Name);
    [~,node_t_id] = ismember(end_nodes(:,2), G.Nodes.Name);

    for outage = 1:width(lodf);
        i = node_s_id(outage);
        j = node_t_id(outage);

        szamlalo = n2sPTDF(:,i) - n2sPTDF(:,j);
        nevezo = 1.0 - szamlalo(outage);
        
        lodf(:,outage) = szamlalo/nevezo;
        lodf(outage,outage) = -1;
    end

end