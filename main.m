% main kód, futtatáshoz, többi függvény meghívásához
% csomópontok és vezetékek tábláinak definiálása
% digraph


clc
clearvars


% Csomópontok
node_name = ["N1"; "N2"; "N3"; "N4"; "N5"; "N6"; "N7"; 
    "N8"; "N9"; "N10"; "N11"; "N12"; "N13"; "N14"; "N15"];

node_ref_export = [114.6666; 55.2255;-61.7903; -42.5794; 156.4802; -61.8938; -95.6438; 
    71.0394;  83.0844; -45.1048; -114.9521; 56.8103; -115.3422; 112.0000; -112.0000];

node_zone = ["A";"A";"B";"B";"A";"B";"B";"C";"C";"C";"C";"D";"D";"D";"D"];

Node_table = table(node_name, node_ref_export, node_zone,'VariableNames',{'Name', 'Reference export', 'Zone'})


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
G = digraph(Edge_table,Node_table);
% [~,loc] = ismember(Edge_table.Name, G.Edges.Name);
% G.Edges(loc,:);
p = plotgraph(G,G.Nodes.Name, G.Edges.Name);

slack = "N10";
zones = unique(G.Nodes.Zone);
% GSK
gsk = zeros(height(G.Nodes), numel(zones));
for col = 1:width(gsk)
    zone_id = col;
    node_id = find(G.Nodes.Zone == zones(zone_id));
    gsk(node_id,col) = 1 / numel(node_id);
end
gsk(find(G.Nodes.Name == slack),:) = [];


% n2sPTDF, z2sPTDF függvény hívása
[n2sPTDF, z2sPTDF] = ptdf(G, Edge_table, slack, gsk);


% Fref
ref_exp = G.Nodes.("Reference export");
ref_exp(find(G.Nodes.Name == slack)) = [];
f_ref = n2sPTDF * ref_exp;
[~,loc]=ismember(G.Edges.Name, Edge_table.Name);
G.Edges.("Reference flow") = f_ref(loc);

p_fref = plotgraph(G, G.Nodes.("Reference export"), G.Edges.("Reference flow"))

% lodf függvény hívása
LODF = lodf(n2sPTDF, slack, G, Edge_table);

% zero balance nettó pozíciók (zónák)
ref_net_pos = zeros(numel(zones),1);
for id = 1:height(ref_net_pos);
    node_id = find(G.Nodes.Zone == zones(id));
    ref_net_pos(id) = sum(G.Nodes.("Reference export")(node_id));
end

% zero balance áramlás
f_0 = f_ref - z2sPTDF * ref_net_pos;
[~,loc]=ismember(G.Edges.Name, Edge_table.Name);
G.Edges.("Zero balance reference flow") = f_0(loc);

% zero balance áramlás plot
% zero balance állapotbeli exportok: ref_exp_0
ref_exp_0 = ref_exp - gsk * ref_net_pos;
slack_id = find(G.Nodes.Name == slack);
ref_exp_0 = [ref_exp_0(1:slack_id-1); 0; ref_exp_0(slack_id:end)];
G.Nodes.("Zero balance reference export") = ref_exp_0;

p_f0 = plotgraph(G, G.Nodes.("Zero balance reference export"), G.Edges.("Zero balance reference flow"))

% f_max, min_ram, frm
f_max = 200;
min_ram = 0.7 * f_max;
frm = 0.1 * f_max;

% cbco előállítása
CBCO_table = cbco(z2sPTDF, LODF, G, Edge_table, f_max, min_ram, frm, zones, ref_net_pos);

% presolve cbco
presolved_CBCO_table = presolve(CBCO_table, zones);