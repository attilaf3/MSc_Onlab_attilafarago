function [n2sPTDF, z2sPTDF] = ptdf(G, Edge_table, slack, gsk)

    % Illeszkedési mátrix
    % full-al a sparse-double-t mátrixba rendezzük
    % - szorzó, mert a matlab alapértelmezetten fordítva értelmezi (a
    % source-hoz rendel -1-et)
    A = incidence(G);
    [~, loc] = ismember(Edge_table.Name, G.Edges.Name);
    A = -full(A(:,loc));
    A = A';
    % slack eliminálása
    A(:,find(G.Nodes.Name == slack)) = [];
    
    % vezeték admittancia mátrixa:
    [~,loc]=ismember(Edge_table.Name, G.Edges.Name);
    Y_line = diag(1./G.Edges.Resistance(loc,:));
    
    % rendszer mátrixa:
    Y_bus = A' * Y_line * A;
    
    % n2sPTDF
    n2sPTDF = Y_line * A / Y_bus;
    
    % z2sPTDF
    z2sPTDF = n2sPTDF * gsk;

end