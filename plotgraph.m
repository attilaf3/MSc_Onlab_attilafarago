function p_ = plotgraph(graph, nodelabel, edgelabel)

    figure
    p_ = plot(graph,'NodeLabel',nodelabel, 'EdgeLabel',edgelabel)
    p_.XData = 2* [0.0; 0.0; 0.0; 0.0; 2.6; 3.4; 3.4; 4.7; 
        5.1; 7.5; 9.7; 9.7; 6.3; 6.2; 9.4];
    p_.YData = 2* [7.5; 5.2; 2.7; 0.4; 5.1; 3.0; 0.6; 7.6; 
        5.0; 5.5; 7.2; 3.3; 2.6; 0.6; 0.7];
    p_.MarkerSize = 7;

    % Különböző zónák csomópontjaihoz eltérő színek rendelése
    highlight(p_, graph.Nodes.Zone =="A",'NodeColor','red');
    highlight(p_, graph.Nodes.Zone =="B",'NodeColor','green');
    highlight(p_, graph.Nodes.Zone =="C",'NodeColor','blue');
    highlight(p_, graph.Nodes.Zone =="D",'NodeColor','magenta');

end