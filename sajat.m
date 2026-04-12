clc
clear all


% Csomóponti table
nodesname = ["N1"; "N2"; "N3"; "N4"; "N5"; "N6"; "N7"; "N8"; "N9"; "N10"; "N11"; "N12"; "N13"; "N14"; "N15"];
nodesnp = [190; 70; -80; -60; 200; -40; -130; 0; 120; 130; -150; 60; -170; 40; -120];
nodesnp(end) = -(sum(nodesnp) - nodesnp(end));
nodesx  = 2* [0.0; 0.0; 0.0; 0.0; 2.2; 3.0; 3.0; 4.2; 4.7; 6.8; 7.8; 8.2; 5.8; 6.0; 8.0];
nodesy  = 2* [5.0; 3.7; 2.2; 0.4; 3.4; 2.0; 0.6; 5.1; 3.5; 4.0; 4.7; 2.8; 2.2; 1.0; 0.7];
nodezone = ["A";"A";"B";"B";"A";"B";"B";"C";"C";"C";"C";"D";"D";"D";"D"];

nodes = table(nodesname, nodesnp, nodesx, nodesy, nodezone, 'VariableNames', {'NodeName', 'NP', 'X', 'Y', 'Zone'});

% Vezeték table
edgesname = ["L1"; "L2"; "L3"; "L4"; "L5"; "L6"; "L7"; "L8"; "L9"; "L10"; "L11"; "L12"; "L13"; "L14"; "L15"; "L16"; "L17"; "L18"; "L19"; "L20"; "L21"; "L22"];
edgeslength = [150; 125; 175; 160; 100; 170; 140; 115; 145; 210; 180; 215; 135; 130; 230; 155; 150; 160; 165; 155; 145; 235];
edgesohm = 0.4 * edgeslength;
edgesfrom = ["N1"; "N2"; "N3"; "N4"; "N6"; "N3"; "N3"; "N2"; "N5"; "N1"; "N8"; "N8"; "N9"; "N10"; "N11"; "N10"; "N9"; "N13"; "N6"; "N7"; "N14"; "N12"];
edgesto   = ["N2"; "N3"; "N4"; "N7"; "N7"; "N6"; "N5"; "N5"; "N9"; "N8"; "N9"; "N11"; "N10"; "N11"; "N12"; "N12"; "N13"; "N12"; "N13"; "N14"; "N15"; "N15"];

edges = table(edgesname, edgesohm, edgesfrom, edgesto, 'VariableNames', {'EdgeName', 'Ohm', 'From', 'To'});

% Slack csomópont
slacknodename = "N10";
slacknodeid = find(nodes.NodeName == slacknodename);
sumnp = sum(nodes.NP) - nodes{slacknodeid, "NP"};
slacknodenp = -sumnp;
nodes.NP(slacknodeid) = slacknodenp;

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
A_table = array2table(A, "RowNames", edgesname, "VariableNames", noslacknodesname)

% Vezeték admittancia mátrix
Yline = diag(1./edges.Ohm);
Yline_table = array2table(Yline, "RowNames", edgesname, "VariableNames", edgesname)

% Rendszer mátrixa mátrixszorzással (slack eliminálva)
Y = A'*Yline*A;
Y_table = array2table(Y, "VariableNames", noslacknodesname, "RowNames", noslacknodesname)

% Node-to-Slack PTDF
nodetoslackPTDF = Yline * A / Y;
nodetoslackPTDF_table = array2table(nodetoslackPTDF, 'VariableNames', noslacknodesname, 'RowNames', edgesname)

% Referencia flow (a nem slack csomópontokra)
np = nodes.NP;
np(slacknodeid) = [];
Fr = nodetoslackPTDF * np;
edges.Flowref = Fr;

% Node-to-Node PTDF (a nem slack csomópontokra)
nodefrom = "N1";
nodeto = "N5";

fromid = find(noslacknodesname == nodefrom);
toid   = find(noslacknodesname == nodeto);

ntncolname = sprintf("%s->%s", nodefrom, nodeto);
nodetonodePTDF = nodetoslackPTDF(:, fromid) - nodetoslackPTDF(:, toid);
nodetonodePTDF_table = array2table(nodetonodePTDF, "RowNames", edgesname, "VariableNames", ntncolname)

% GSK (egyenletes)
zonenames = ["A","B","C","D"];
zonenum = numel(zonenames);

GSK = zeros(height(nodes), zonenum);

for k = 1:zonenum
    nodeid = find(nodes.Zone == zonenames(k));
    GSK(nodeid,k) = 1 / numel(nodeid);
end

GSK(slacknodeid,:) = [];
GSK_table = array2table(GSK, "RowNames", noslacknodesname,"VariableNames", zonenames)

% Zone-to-Slack PTDF
zoneplot = "B";
zoneid_plot = find(zonenames == zoneplot);

zonetoslackPTDF = nodetoslackPTDF * GSK;
zonetoslackPTDF_table = array2table(zonetoslackPTDF, "RowNames", edgesname, "VariableNames", zonenames)

% Zone-to-Zone PTDF
zonefrom = "A";
zoneto   = "B";
fromid = find(zonenames == zonefrom);
toid   = find(zonenames == zoneto);

zonetozonePTDF = nodetoslackPTDF * (GSK(:,fromid) - GSK(:,toid));
ztzcolname = sprintf("%s->%s", zonefrom, zoneto);
ztzPTDF_table = array2table(zonetozonePTDF,"RowNames", edgesname,"VariableNames", ztzcolname)


% LODF
outageline = "L22";
outageid = find(edges.EdgeName == outageline);
LODF = zeros(height(edges),height(edges));
nodetoslackPTDF_full = [nodetoslackPTDF(:,1:slacknodeid-1), zeros(size(nodetoslackPTDF,1), 1), nodetoslackPTDF(:,slacknodeid:end)];
for k = 1:height(edges)

    from = edges.From(k);
    to   = edges.To(k);

    fromid = find(nodesname == from);
    toid   = find(nodesname == to);

    ptdf_ij = nodetoslackPTDF_full(:,fromid) - nodetoslackPTDF_full(:,toid);
    LODF(:,k) = ptdf_ij / (1 - ptdf_ij(k));
    LODF(k,k) = -1;
end


LODF_table = array2table(LODF, "RowNames", edgesname, "VariableNames", edgesname)
% oszlopban az adott kieső vezeték
% sorban a kieső vezeték hatása az adott vezetékre


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

function plot_network(nodes, edges, slacknodename, figTitle)

    G = digraph(cellstr(edges.From), cellstr(edges.To), ones(height(edges),1), cellstr(nodes.NodeName));

    figure('Color','w')
    p = plot(G, ...
        'XData', nodes.X, ...
        'YData', nodes.Y, ...
        'LineWidth', 1.5, ...
        'ArrowSize', 12);        

    p.NodeLabel = repmat({''}, height(nodes),1);
    p.EdgeLabel = cellstr(edges.EdgeName);
    
    axis equal
    axis off
    title(figTitle, 'Interpreter','none')

    hold on

    nodeColors = get_zone_colors(nodes.Zone);

    for node = 1:height(nodes)
        if nodes.NodeName(node) == slacknodename
            edgeColor = [0 0 0];
            lw = 2.2;
        else
            edgeColor = [0.35 0.35 0.35];
            lw = 1.2;
        end

        scatter(nodes.X(node), nodes.Y(node), 2300, ...
            'MarkerFaceColor', nodeColors(node,:), ...
            'MarkerEdgeColor', edgeColor, ...
            'LineWidth', lw);

        txt = sprintf('%s\nNP=%+.1f', nodes.NodeName(node), nodes.NP(node));
        text(nodes.X(node), nodes.Y(node)+0.10, txt, ...
            'HorizontalAlignment','center', ...
            'FontWeight','bold', ...
            'FontSize',9, ...
            'Interpreter','none');
    end
    
    hold off
end

function plot_edge_values(nodes, edges, edgeValues, slacknodename, figTitle)

    G = digraph(cellstr(edges.From), cellstr(edges.To), ones(height(edges),1), cellstr(nodes.NodeName));

    figure('Color','w')
    p = plot(G, ...
        'XData', nodes.X, ...
        'YData', nodes.Y, ...
        'LineWidth', 1.5, ...
        'ArrowSize', 12);

    p.NodeLabel = repmat({''}, height(nodes),1);
    edgeLabels = compose("%s  %+.4f", string(edges.EdgeName), edgeValues);
    p.EdgeLabel = cellstr(edgeLabels);

    axis equal
    axis off
    title(figTitle, 'Interpreter','none')

    hold on

    nodeColors = get_zone_colors(nodes.Zone);

    for node = 1:height(nodes)
        if nodes.NodeName(node) == slacknodename
            edgeColor = [0 0 0];
            lw = 2.2;
        else
            edgeColor = [0.35 0.35 0.35];
            lw = 1.2;
        end

        scatter(nodes.X(node), nodes.Y(node), 2300, ...
            'MarkerFaceColor', nodeColors(node,:), ...
            'MarkerEdgeColor', edgeColor, ...
            'LineWidth', lw);

        txtNode = sprintf('%s\nNP=%+.1f', nodes.NodeName(node), nodes.NP(node));
        text(nodes.X(node), nodes.Y(node)+0.10, txtNode, ...
            'HorizontalAlignment','center', ...
            'FontWeight','bold', ...
            'FontSize',9, ...
            'Interpreter','none');
    end
    
    hold off
end

function plot_lodf_case(nodes, edges, lodfVec, outageid, slacknodename, figTitle)

    G = digraph(cellstr(edges.From), cellstr(edges.To), ones(height(edges),1), cellstr(nodes.NodeName));

    figure('Color','w')
    p = plot(G, ...
        'XData', nodes.X, ...
        'YData', nodes.Y, ...
        'LineWidth', 1.5, ...
        'ArrowSize', 12);

    p.NodeLabel = repmat({''}, height(nodes), 1);
    edgeLabels = compose("%s  %+.4f", string(edges.EdgeName), lodfVec);
    edgeLabels(outageid) = string(edges.EdgeName(outageid)) + " KIESIK";
    p.EdgeLabel = cellstr(edgeLabels);

    axis equal
    axis off
    title(figTitle, 'Interpreter','none')

    hold on

    nodeColors = get_zone_colors(nodes.Zone);

    for node = 1:height(nodes)
        if nodes.NodeName(node) == slacknodename
            edgeColor = [0 0 0];
            lw = 2.2;
        else
            edgeColor = [0.35 0.35 0.35];
            lw = 1.2;
        end

        scatter(nodes.X(node), nodes.Y(node), 2300, ...
            'MarkerFaceColor', nodeColors(node,:), ...
            'MarkerEdgeColor', edgeColor, ...
            'LineWidth', lw);

        txtNode = sprintf('%s\nNP=%+.1f', char(nodes.NodeName(node)), nodes.NP(node));
        text(nodes.X(node), nodes.Y(node)+0.10, txtNode, ...
            'HorizontalAlignment','center', ...
            'FontWeight','bold', ...
            'FontSize',9, ...
            'Interpreter','none');
    end

    hold off
end

function nodeColors = get_zone_colors(zones)

    nodeColors = zeros(numel(zones), 3);

    for k = 1:numel(zones)
        switch char(zones(k))
            case 'A'
                nodeColors(k,:) = [0.85 0.20 0.20];
            case 'B'
                nodeColors(k,:) = [0.15 0.60 0.25];
            case 'C'
                nodeColors(k,:) = [0.50 0.20 0.70];
            case 'D'
                nodeColors(k,:) = [0.95 0.35 0.55];
            otherwise
                nodeColors(k,:) = [0.5 0.5 0.5];
        end
    end
end

%[appendix]
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
