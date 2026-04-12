function results = matrices(nodes, edges, slacknodeid)

noslacknodesname = nodes.NodeName;
noslacknodesname(slacknodeid) = [];

% Illeszkedési mátrix
A = zeros(height(edges), height(nodes));

for k = 1:height(edges)
    fromid = find(nodes.NodeName == edges.From(k));
    toid = find(nodes.NodeName == edges.To(k));
    A(k, fromid) = 1;
    A(k, toid) = -1;
end

A(:,slacknodeid) = [];
A_table = array2table(A, "RowNames", edges.EdgeName, "VariableNames", noslacknodesname);

% Vezeték admittancia mátrix
Yline = diag(1./edges.Ohm);
Yline_table = array2table(Yline, "RowNames", edges.EdgeName, "VariableNames", edges.EdgeName);

% Rendszer mátrixa mátrixszorzással (slack eliminálva)
Y = A' * Yline * A;
Y_table = array2table(Y, "VariableNames", noslacknodesname, "RowNames", noslacknodesname);

% Node-to-Slack PTDF
nodetoslackPTDF = Yline * A / Y;
nodetoslackPTDF_table = array2table(nodetoslackPTDF, ...
    'VariableNames', noslacknodesname, 'RowNames', edges.EdgeName);

% Referencia flow (a nem slack csomópontokra)
np = nodes.NP;
np(slacknodeid) = [];
Fr = nodetoslackPTDF * np;

% Node-to-Node PTDF (a nem slack csomópontokra)
nodefrom = "N1";
nodeto = "N5";

fromid = find(noslacknodesname == nodefrom);
toid   = find(noslacknodesname == nodeto);

ntncolname = sprintf("%s->%s", nodefrom, nodeto);
nodetonodePTDF = nodetoslackPTDF(:, fromid) - nodetoslackPTDF(:, toid);
nodetonodePTDF_table = array2table(nodetonodePTDF, ...
    "RowNames", edges.EdgeName, "VariableNames", ntncolname);

% GSK (egyenletes)
zonenames = ["A","B","C","D"];
zonenum = numel(zonenames);

GSK = zeros(height(nodes), zonenum);

for k = 1:zonenum
    nodeid = find(nodes.Zone == zonenames(k));
    GSK(nodeid,k) = 1 / numel(nodeid);
end

GSK(slacknodeid,:) = [];
GSK_table = array2table(GSK, "RowNames", noslacknodesname, "VariableNames", zonenames);

% Zone-to-Slack PTDF
zoneplot = "B";
zoneid_plot = find(zonenames == zoneplot);

zonetoslackPTDF = nodetoslackPTDF * GSK;
zonetoslackPTDF_table = array2table(zonetoslackPTDF, ...
    "RowNames", edges.EdgeName, "VariableNames", zonenames);

% Zone-to-Zone PTDF
zonefrom = "A";
zoneto   = "B";
fromid = find(zonenames == zonefrom);
toid   = find(zonenames == zoneto);

zonetozonePTDF = nodetoslackPTDF * (GSK(:,fromid) - GSK(:,toid));
ztzcolname = sprintf("%s->%s", zonefrom, zoneto);
ztzPTDF_table = array2table(zonetozonePTDF, ...
    "RowNames", edges.EdgeName, "VariableNames", ztzcolname);

% LODF
outageline = "L22";
outageid = find(edges.EdgeName == outageline);
LODF = zeros(height(edges), height(edges));
nodetoslackPTDF_full = [nodetoslackPTDF(:,1:slacknodeid-1), ...
                        zeros(size(nodetoslackPTDF,1), 1), ...
                        nodetoslackPTDF(:,slacknodeid:end)];

for k = 1:height(edges)

    from = edges.From(k);
    to   = edges.To(k);

    fromid = find(nodes.NodeName == from);
    toid   = find(nodes.NodeName == to);

    ptdf_ij = nodetoslackPTDF_full(:,fromid) - nodetoslackPTDF_full(:,toid);
    LODF(:,k) = ptdf_ij / (1 - ptdf_ij(k));
    LODF(k,k) = -1;
end

LODF_table = array2table(LODF, "RowNames", edges.EdgeName, "VariableNames", edges.EdgeName);

% Referencia flow tárolása külön edge táblában
edges_with_flow = edges;
edges_with_flow.Flowref = Fr;

% Kimenet
results = struct();

results.A = A;
results.A_table = A_table;

results.Yline = Yline;
results.Yline_table = Yline_table;

results.Y = Y;
results.Y_table = Y_table;

results.nodetoslackPTDF = nodetoslackPTDF;
results.nodetoslackPTDF_table = nodetoslackPTDF_table;

results.Fr = Fr;
results.edges_with_flow = edges_with_flow;

results.nodefrom = nodefrom;
results.nodeto = nodeto;
results.nodetonodePTDF = nodetonodePTDF;
results.nodetonodePTDF_table = nodetonodePTDF_table;

results.zonenames = zonenames;
results.GSK = GSK;
results.GSK_table = GSK_table;

results.zoneplot = zoneplot;
results.zoneid_plot = zoneid_plot;
results.zonetoslackPTDF = zonetoslackPTDF;
results.zonetoslackPTDF_table = zonetoslackPTDF_table;

results.zonefrom = zonefrom;
results.zoneto = zoneto;
results.zonetozonePTDF = zonetozonePTDF;
results.ztzPTDF_table = ztzPTDF_table;

results.outageline = outageline;
results.outageid = outageid;
results.nodetoslackPTDF_full = nodetoslackPTDF_full;
results.LODF = LODF;
results.LODF_table = LODF_table;

end

%[appendix]
%---
