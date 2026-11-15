% main kód, futtatáshoz, többi függvény meghívásához
% csomópontok és vezetékek tábláinak definiálása
% digraph


clc
clearvars


% Csomópontok
node_name = ["N1"; "N2"; "N3"; "N4"; "N5"; "N6"; "N7"; 
    "N8"; "N9"; "N10"; "N11"; "N12"; "N13"; "N14"; "N15"];

node_ref_balance = [114.6666; 55.2255;-61.7903; -42.5794; 156.4802; -61.8938; -95.6438; 
    71.0394;  83.0844; -45.1048; -114.9521; 56.8103; -115.3422; 112.0000; -112.0000];

node_zone = ["A";"A";"B";"B";"A";"B";"B";"C";"C";"C";"C";"D";"D";"D";"D"];

Node_table = table(node_name, node_ref_balance, node_zone,'VariableNames',{'Name', 'Reference balance', 'Zone'})


% Vezetékek
edge_name = ["L1"; "L2"; "L3"; "L4"; "L5"; "L6"; "L7"; "L8"; "L9"; "L10"; "L11";
    "L12"; "L13"; "L14"; "L15"; "L16"; "L17"; "L18"; "L19"; "L20"; "L21"; "L22"];

edge_km = [150; 125; 175; 160; 100; 170; 140; 115; 145; 210; 180; 
    215; 135; 130; 230; 155; 150; 160; 165; 155; 145; 235];

edge_ohm = 0.3 * edge_km;

edge_s = ["N1"; "N2"; "N3"; "N4"; "N6"; "N3"; "N3"; "N2"; "N5"; "N1"; "N8";
    "N8"; "N9"; "N10"; "N11"; "N10"; "N9"; "N13"; "N6"; "N7"; "N14"; "N12"];

edge_t = ["N2"; "N3"; "N4"; "N7"; "N7"; "N6"; "N5"; "N5"; "N9"; "N8"; "N9";
    "N11"; "N10"; "N11"; "N12"; "N12"; "N13"; "N12"; "N13"; "N14"; "N15"; "N15"];

Edge_table = table([edge_s edge_t], edge_name, edge_ohm, ...
    'VariableNames',{'EndNodes', 'Name','Resistance'})


% Digraph
G = digraph(Edge_table,Node_table)

p = plot(G,'NodeLabel',G.Nodes.Name, 'EdgeLabel',G.Edges.Name)
p.XData = 2* [0.0; 0.0; 0.0; 0.0; 2.6; 3.4; 3.4; 4.7; 
    5.1; 7.5; 9.7; 9.7; 6.3; 6.2; 9.4];
p.YData = 2* [7.5; 5.2; 2.7; 0.4; 5.1; 3.0; 0.6; 7.6; 
    5.0; 5.5; 7.2; 3.3; 2.6; 0.6; 0.7];
p.MarkerSize = 7;

% Különböző zónák csomópontjaihoz eltérő színek rendelése
highlight(p, G.Nodes.Zone =="A",'NodeColor','red')
highlight(p, G.Nodes.Zone =="B",'NodeColor','green')
highlight(p, G.Nodes.Zone =="C",'NodeColor','blue')
highlight(p, G.Nodes.Zone =="D",'NodeColor','magenta')



% ez külön függvény: n2s_PTDF
% Illeszkedési mátrix
% full-al a sparse-double-t mátrixba rendezzük
% - szorzó, mert a matlab alapértelmezetten fordítva értelmezi (a
% source-hoz rendel -1-et)
A = incidence(G);
[~, loc] = ismember(Edge_table.Name, G.Edges.Name);
A = full(A(:,loc));
% slack eliminálása
slack = "N10";
A(find(G.Nodes.Name == slack),:) = [];

% vezeték admittancia mátrixa:
[~,loc]=ismember(Edge_table.Name, G.Edges.Name);
Y_line = diag(1./G.Edges.Resistance(loc,:));

% rendszer mátrixa:
Y_bus = A * Y_line * A';


