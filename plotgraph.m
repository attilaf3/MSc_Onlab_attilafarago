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
