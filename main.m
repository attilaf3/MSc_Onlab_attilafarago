clc
clear all

% Hálózat összeállítása
[nodes, edges, slacknodename, slacknodeid, noslacknodesname] = build_network_data();

% Mátrixok és számolt mennyiségek
results = matrices(nodes, edges, slacknodeid);


A_table = results.A_table;
Yline_table = results.Yline_table;
Y_table = results.Y_table;
nodetoslackPTDF_table = results.nodetoslackPTDF_table;
nodetonodePTDF_table = results.nodetonodePTDF_table;
GSK_table = results.GSK_table;
zonetoslackPTDF_table = results.zonetoslackPTDF_table;
ztzPTDF_table = results.ztzPTDF_table;
LODF_table = results.LODF_table;

Fr = results.Fr;
edges = results.edges_with_flow;

nodefrom = results.nodefrom;
nodeto = results.nodeto;
nodetonodePTDF = results.nodetonodePTDF;

zoneplot = results.zoneplot;
zoneid_plot = results.zoneid_plot;
zonetoslackPTDF = results.zonetoslackPTDF;

zonefrom = results.zonefrom;
zoneto = results.zoneto;
zonetozonePTDF = results.zonetozonePTDF;

outageline = results.outageline;
outageid = results.outageid;
LODF = results.LODF;

% Hálózat alaprajza
plot_network(nodes, edges, slacknodename, 'Mintahálózat');

% Eredeti (referencia) áramlás az éleken
plot_edge_values(nodes, edges, Fr, slacknodename, ...
    'Referencia áramlás az éleken');

% Node-to-node PTDF
mw_ntn = 1;
flow_ntn = mw_ntn * nodetonodePTDF;
plot_edge_values(nodes, edges, flow_ntn, slacknodename, ...
    sprintf('Áramlás + %g MW node-to-node: %s -> %s', ...
    mw_ntn, nodefrom, nodeto));

% Node-to-node PTDF után flow
mw_ntn = 1;
flow_ntn = Fr + mw_ntn * nodetonodePTDF;
plot_edge_values(nodes, edges, flow_ntn, slacknodename, ...
    sprintf('Áramlás + %g MW node-to-node: %s -> %s', ...
    mw_ntn, nodefrom, nodeto));

% Zone-to-slack PTDF egy kiválasztott zónára
plot_edge_values(nodes, edges, zonetoslackPTDF(:, zoneid_plot), slacknodename, ...
    sprintf('Zone-to-slack PTDF: %s -> slack (1 MW)', zoneplot));

% Zone-to-slack áramlás egy kiválasztott zónára
mw_zts = 1;
flow_z2s = Fr + mw_zts * zonetoslackPTDF(:, zoneid_plot);
plot_edge_values(nodes, edges, flow_z2s, slacknodename, ...
    sprintf('Áramlás + %g MW zone-to-slack: %s -> slack', ...
    mw_zts, zoneplot));

% Zone-to-zone PTDF adott zónából egy másik, megadott zónába
plot_edge_values(nodes, edges, zonetozonePTDF, slacknodename, ...
    sprintf('Zone-to-zone PTDF: %s -> %s (1 MW)', zonefrom, zoneto));

% Zone-to-zone áramlás adott zónából egy másik, megadott zónába
mw_ztz = 1;
flow_z2z = Fr + mw_ztz * zonetozonePTDF;
plot_edge_values(nodes, edges, flow_z2z, slacknodename, ...
    sprintf('Áramlás + %g MW zone-to-zone: %s -> %s', ...
    mw_ztz, zonefrom, zoneto));

% LODF egy kiválasztott kieső vezetékre
plot_lodf_case(nodes, edges, LODF(:, outageid), outageid, slacknodename, ...
    sprintf('LODF együtthatók %s kiesése esetén', outageline));

% Kiesés utáni flow
flowaftout = edges.Flowref + LODF(:, outageid) * edges.Flowref(outageid);
flowaftout(outageid) = 0;

plot_edge_values(nodes, edges, flowaftout, slacknodename, ...
    sprintf('%s kiesése utáni áramlások', outageline));


%[appendix]
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
