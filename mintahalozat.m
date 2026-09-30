clc
clear all

%% Configure network
% Csomóponti table
nodesname = ["N1"; "N2"; "N3"; "N4"; "N5"; "N6"; "N7"; "N8"; "N9"; "N10"; "N11"; "N12"; "N13"; "N14"; "N15"];
nodesnp = [190; 70; -80; -60; 200; -90; -130; 110; 120; 130; -150; 60; -170; 170; -120];
nodesnp(end) = -(sum(nodesnp) - nodesnp(end));
nodesx  = 2* [0.0; 0.0; 0.0; 0.0; 2.6; 3.4; 3.4; 4.7; 5.1; 7.5; 9.7; 9.7; 6.3; 6.2; 9.4];
nodesy  = 2* [7.5; 5.2; 2.7; 0.4; 5.1; 3.0; 0.6; 7.6; 5.0; 5.5; 7.2; 3.3; 2.6; 0.6; 0.7];
nodezone = ["A";"A";"B";"B";"A";"B";"B";"C";"C";"C";"C";"D";"D";"D";"D"];


% lsqlin-ből kapott np, ha finalram
% nodesnp = 0.8 * [ ...
%      143.3333;
%       69.0319;
%      -77.2379;
%      -53.2243;
%      195.6002;
%      -77.3673;
%     -119.5547;
%       88.7992;
%      103.8555;
%      -56.3812;
%     -143.6901;
%       71.0129;
%     -144.1777;
%      140.0000;
%     -140.0000
% ];
nodesnp = [...
  114.6666;
   55.2255;
  -61.7903;
  -42.5794;
  156.4802;
  -61.8938;
  -95.6438;
   71.0394;
   83.0844;
  -45.1048;
 -114.9521;
   56.8103;
 -115.3422;
  112.0000;
 -112.0000
    ];



nodes = table(nodesname, nodesnp, nodesx, nodesy, nodezone, 'VariableNames', {'NodeName', 'NP', 'X', 'Y', 'Zone'});

% Vezeték table
edgesname = ["L1"; "L2"; "L3"; "L4"; "L5"; "L6"; "L7"; "L8"; "L9"; "L10"; "L11"; "L12"; "L13"; "L14"; "L15"; "L16"; "L17"; "L18"; "L19"; "L20"; "L21"; "L22"];
edgeslength = [150; 125; 175; 160; 100; 170; 140; 115; 145; 210; 180; 215; 135; 130; 230; 155; 150; 160; 165; 155; 145; 235];
edgesohm = 0.3 * edgeslength;
edgesfrom = ["N1"; "N2"; "N3"; "N4"; "N6"; "N3"; "N3"; "N2"; "N5"; "N1"; "N8"; "N8"; "N9"; "N10"; "N11"; "N10"; "N9"; "N13"; "N6"; "N7"; "N14"; "N12"];
edgesto   = ["N2"; "N3"; "N4"; "N7"; "N7"; "N6"; "N5"; "N5"; "N9"; "N8"; "N9"; "N11"; "N10"; "N11"; "N12"; "N12"; "N13"; "N12"; "N13"; "N14"; "N15"; "N15"];
% edgesfmax = 1386 * ones(numel(edgesname),1);
edgesfmax = 200 * ones(numel(edgesname),1);
edgesfrm  = 0.1*edgesfmax;   

edges = table(edgesname, edgesohm, edgesfrom, edgesto, edgesfmax, edgesfrm, 'VariableNames', {'EdgeName', 'Ohm', 'From', 'To', 'Fmax', 'FRM'});

% Slack csomópont
slacknodename = "N10";
slacknodeid = find(nodes.NodeName == slacknodename);

noslacknodesname = nodes.NodeName;
noslacknodesname(slacknodeid) = [];

%% Build A, Y, PTDF matrices
% Illeszkedési mátrix
A = zeros(height(edges), height(nodes));

for edge_idx = 1:height(edges)
    src_nodeidx = find(nodes.NodeName == edges.From(edge_idx));
    dst_nodeidx = find(nodes.NodeName == edges.To(edge_idx));
    A(edge_idx, src_nodeidx) = 1;
    A(edge_idx, dst_nodeidx) = -1;
end


% Slacknode oszlopának törlése (DC loadflowban a fázisszög 0-nak való
% rögzítése)
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

% Node-to-Node PTDF (a nem slack csomópontokra)
nodefrom = "N1";
nodeto = "N5";

% fromid = find(noslacknodesname == nodefrom);
% toid   = find(noslacknodesname == nodeto);

fromid = noslacknodesname == nodefrom;
toid   = noslacknodesname == nodeto;

ntncolname = sprintf("%s->%s", nodefrom, nodeto);
nodetonodePTDF = nodetoslackPTDF(:, fromid) - nodetoslackPTDF(:, toid);
nodetonodePTDF_table = array2table(nodetonodePTDF, "RowNames", edgesname, "VariableNames", ntncolname)

%% Create reference flow
np = nodes.NP;
np(slacknodeid) = [];
edges.Fref  = nodetoslackPTDF * np;

%% Create GSK (egyenletes) and zone PTDFs
zonenames = ["A","B","C","D"];
zonenum = numel(zonenames);

GSK = zeros(height(nodes), zonenum);

for edge_idx = 1:zonenum
    nodeid = find(nodes.Zone == zonenames(edge_idx));
    GSK(nodeid,edge_idx) = 1 / numel(nodeid);
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

% Zone-to-Zone PTDF minden kombinációra

ztz_pairs = [];
ztz_names = strings(0,1);

k = 0;
for a = 1:zonenum
    for b = 1:zonenum
        if a ~= b
            k = k + 1;
            ztz_pairs(k,:) = [a b];
            ztz_names(k) = zonenames(a) + "->" + zonenames(b);
        end
    end
end

ztznum = size(ztz_pairs,1);
zonetozonePTDF_all = zeros(height(edges), ztznum);

for k = 1:ztznum
    fromZ = ztz_pairs(k,1);
    toZ   = ztz_pairs(k,2);

    zonetozonePTDF_all(:,k) = ...
        zonetoslackPTDF(:,fromZ) - zonetoslackPTDF(:,toZ);
end

zonetozonePTDF_table = array2table( ...
    zonetozonePTDF_all, ...
    "RowNames", edgesname, ...
    "VariableNames", ztz_names);


% Egyes zónákban a nettó import/export (zónás nettó pozíció)
zoneNP = zeros(zonenum,1);

for edge_idx = 1:zonenum
    nodeid = find(nodes.Zone == zonenames(edge_idx));
    zoneNP(edge_idx) = sum(nodes.NP(nodeid));
end

zoneNP_table = table(zonenames', zoneNP, ...
    'VariableNames', {'Zone', 'ZoneNP_ref'});

%% ZERO BALANCE és VIRTUÁLIS KAPACITÁS SZÁMÍTÁS
% Kereskedelem nélküli maradékáramlás:
% F0 = Fref - (nodetoslackPTDF * GSK) * zoneNP
% F0 = Fref - zonetoslackPTDF * zoneNP
edges.F0 = edges.Fref - zonetoslackPTDF * zoneNP;

% RAM
edges.RAM0 = edges.Fmax - edges.FRM - edges.F0;


% Minimum RAM követelmény
Ramr = 0.70;
edges.minRAM = Ramr * edges.Fmax;

% Pozitív irány
% PTDF*NP <= RAM_plus
edges.RAM0_plus = edges.Fmax - edges.FRM - edges.F0;
edges.AMR_plus = max(edges.minRAM - edges.RAM0_plus, 0);
edges.finalRAM_plus = edges.RAM0_plus + edges.AMR_plus;

% Negatív irány
% -PTDF*NP <= RAM_minus
edges.RAM0_minus = edges.Fmax - edges.FRM + edges.F0;
edges.AMR_minus = max(edges.minRAM - edges.RAM0_minus, 0);
edges.finalRAM_minus = edges.RAM0_minus + edges.AMR_minus;

ram_table = edges(:, {'EdgeName','Fref','F0','Fmax','FRM', ...
    'RAM0_plus','AMR_plus','finalRAM_plus', ...
    'RAM0_minus','AMR_minus','finalRAM_minus','minRAM'})

%% LODF
outageline = "L22";
outageid = find(edges.EdgeName == outageline);
LODF = zeros(height(edges),height(edges));
nodetoslackPTDF_full = [nodetoslackPTDF(:,1:slacknodeid-1), zeros(size(nodetoslackPTDF,1), 1), nodetoslackPTDF(:,slacknodeid:end)];
for edge_idx = 1:height(edges)

    srcnode = edges.From(edge_idx);
    dstnode = edges.To(edge_idx);

    fromid = nodesname == srcnode;
    toid   = nodesname == dstnode;

    ptdf_ij = nodetoslackPTDF_full(:,fromid) - nodetoslackPTDF_full(:,toid);
    LODF(:,edge_idx) = ptdf_ij / (1 - ptdf_ij(edge_idx));
    LODF(edge_idx,edge_idx) = -1;
end


LODF_table = array2table(LODF, "RowNames", edgesname, "VariableNames", edgesname)
% oszlopban az adott kieső vezeték
% sorban a kieső vezeték hatása az adott vezetékre


%% CBCO

critical_outages = ["L1", "L2","L7","L20", "L21","L22"]; 

critical_outages_idx = [];
for k = 1:numel(critical_outages)
    idx = find(edges.EdgeName == critical_outages(k));
    if ~isempty(idx)
        critical_outages_idx(end+1) = idx;
    end
end

cb_names = strings(0,1);
co_names = strings(0,1);      
fref_list = [];
f0_list = [];
ptdf_cbco_matrix = [];
minram_list = [];
ram0_plus_list = [];
ram0_minus_list = [];
amr_plus_list = [];
amr_minus_list = [];
finalram_plus_list = [];
finalram_minus_list = [];
fmax_list = [];


% alapeset
for cb_idx = 1:height(edges)
    cb_name = edges.EdgeName(cb_idx);
    
    cb_names(end+1,1) = cb_name;
    co_names(end+1,1) = "BaseCase";
    
    fref_base = edges.Fref(cb_idx);
    f0_base = edges.F0(cb_idx);

    ram0_plus_base  = edges.Fmax(cb_idx) - edges.FRM(cb_idx) - f0_base;
    ram0_minus_base = edges.Fmax(cb_idx) - edges.FRM(cb_idx) + f0_base;

    minram_base = Ramr * edges.Fmax(cb_idx);

    amr_plus_base  = max(minram_base - ram0_plus_base, 0);
    amr_minus_base = max(minram_base - ram0_minus_base, 0);

    finalram_plus_base  = ram0_plus_base  + amr_plus_base;
    finalram_minus_base = ram0_minus_base + amr_minus_base;
    
    
    ram0_plus_list(end+1,1) = ram0_plus_base;
    ram0_minus_list(end+1,1) = ram0_minus_base;
    minram_list(end+1,1) = minram_base;
    amr_plus_list(end+1,1) = amr_plus_base;
    amr_minus_list(end+1,1) = amr_minus_base;
    finalram_plus_list(end+1,1) = finalram_plus_base;
    finalram_minus_list(end+1,1) = finalram_minus_base;
    f0_list(end+1,1) = f0_base;
    fref_list(end+1,1) = fref_base;
    fmax_list(end+1,1) = edges.Fmax(cb_idx);
    
    % ptdf_cbco_matrix(end+1,:) = zonetozonePTDF_all(cb_idx,:);
    ptdf_cbco_matrix(end+1,:) = zonetoslackPTDF(cb_idx,:);
end

% kontingenciák
for m = 1:numel(critical_outages_idx)
    co_idx = critical_outages_idx(m);
    co_name = edges.EdgeName(co_idx);
    
    for cb_idx = 1:height(edges)
        
        if cb_idx == co_idx
            continue; 
        end
        
        cb_name = edges.EdgeName(cb_idx);
        
        
        
        % ptdf_cbo = zonetozonePTDF_all(cb_idx, :) + LODF(cb_idx, co_idx) * zonetozonePTDF_all(co_idx, :);
        ptdf_cbo = zonetoslackPTDF(cb_idx,:) + LODF(cb_idx, co_idx) * zonetoslackPTDF(co_idx,:);

        fref_cbo = edges.Fref(cb_idx) + LODF(cb_idx, co_idx) * edges.Fref(co_idx);
        
        % f0_cbo = edges.F0(cb_idx) + LODF(cb_idx, co_idx) * edges.F0(co_idx);
        f0_cbo = fref_cbo - ptdf_cbo * zoneNP;



        ram0_plus_cbo  = edges.Fmax(cb_idx) - edges.FRM(cb_idx) - f0_cbo;
        ram0_minus_cbo = edges.Fmax(cb_idx) - edges.FRM(cb_idx) + f0_cbo;
        
        minram_cbo = Ramr * edges.Fmax(cb_idx);
        
        amr_plus_cbo  = max(minram_cbo - ram0_plus_cbo, 0);
        amr_minus_cbo = max(minram_cbo - ram0_minus_cbo, 0);
        
        finalram_plus_cbo  = ram0_plus_cbo  + amr_plus_cbo;
        finalram_minus_cbo = ram0_minus_cbo + amr_minus_cbo;
        
          
        cb_names(end+1,1) = cb_name;
        co_names(end+1,1) = co_name + " kiesése";

        fref_list(end+1,1) = fref_cbo;
        f0_list(end+1,1) = f0_cbo;
        ptdf_cbco_matrix(end+1,:) = ptdf_cbo;


        ram0_plus_list(end+1,1) = ram0_plus_cbo;
        ram0_minus_list(end+1,1) = ram0_minus_cbo;
        minram_list(end+1,1) = minram_cbo;
        amr_plus_list(end+1,1) = amr_plus_cbo;
        amr_minus_list(end+1,1) = amr_minus_cbo;
        finalram_plus_list(end+1,1) = finalram_plus_cbo;
        finalram_minus_list(end+1,1) = finalram_minus_cbo;
        fmax_list(end+1,1) = edges.Fmax(cb_idx);
    end
end

CBCO_table = table(cb_names, co_names, fmax_list, fref_list, ...
    finalram_plus_list, finalram_minus_list, amr_plus_list, amr_minus_list, ...
    'VariableNames', {'CriticalBranch', 'CriticalOutage', ...
    'Fmax', 'Fref', 'RAM_plus', 'RAM_minus', 'AMR_plus', 'AMR_minus'});

% ztzPTDF hozzáadás
% for z = 1:ztznum
%     col_name = "PTDF_" + replace(ztz_names(z), "->", "_to_");
%     CBCO_table.(col_name) = ptdf_cbco_matrix(:, z);
% end
% ztsPTDF hozzáadás
for z = 1:zonenum
    col_name = "PTDF_" + zonenames(z);
    CBCO_table.(col_name) = ptdf_cbco_matrix(:, z);
end

%% CBCO-ba node-to-slack PTDF

ntsPTDF_cbco = [];

% Base case
for cb_idx = 1:height(edges)
    ntsPTDF_cbco(end+1,:) = nodetoslackPTDF(cb_idx,:);
end

% Kontingenciák
for m = 1:numel(critical_outages_idx)
    co_idx = critical_outages_idx(m);

    for cb_idx = 1:height(edges)

        if cb_idx == co_idx
            continue;
        end

        ptdf_cbco_node = nodetoslackPTDF(cb_idx,:) + ...
            LODF(cb_idx, co_idx) * nodetoslackPTDF(co_idx,:);

        ntsPTDF_cbco(end+1,:) = ptdf_cbco_node;
    end
end

% for n = 1:numel(noslacknodesname)
%     col_name = "PTDF_" + noslacknodesname(n);
%     CBCO_table.(col_name) = ntzPTDF_cbco(:,n);
% end
disp(CBCO_table)

%% F0 és RAM0 tábla CBCO-kra

F0_RAM0_table = table( ...
    string(CBCO_table.CriticalBranch), ...
    string(CBCO_table.CriticalOutage), ...
    f0_list, ...
    edges.Fmax(1) * ones(size(f0_list)), ...
    edges.FRM(1)  * ones(size(f0_list)), ...
    ram0_plus_list, ...
    ram0_minus_list, ...
    minram_list, ...
    amr_plus_list, ...
    amr_minus_list, ...
    finalram_plus_list, ...
    finalram_minus_list, ...
    'VariableNames', {'CB','CO','F0','Fmax','FRM', ...
    'RAM0_plus','RAM0_minus','minRAM', ...
    'AMR_plus','AMR_minus','finalRAM_plus','finalRAM_minus'});
%% PTDF * NP <= RAM tábla

% PTDF⋅NP<=RAM+
% −PTDF⋅NP<=RAM−
NP = zoneNP;                 
             
left_plus = ptdf_cbco_matrix * zoneNP;
left_minus = -ptdf_cbco_matrix * zoneNP;

% ha pozitívak, akkor teljesülnek a feltételek
dif_plus  = finalram_plus_list  - left_plus;
dif_minus = finalram_minus_list - left_minus;

PTDF_NP_RAM_table = table( ...
    CBCO_table.CriticalBranch, ...
    CBCO_table.CriticalOutage, ...
    left_plus, finalram_plus_list, dif_plus, ...
    left_minus, finalram_minus_list, dif_minus, ...
    'VariableNames', { ...
        'CriticalBranch', 'CriticalOutage', ...
        'PTDFxNP_plus', 'RAM_plus', 'Difference_plus', ...
        'PTDFxNP_minus', 'RAM_minus', 'Difference_minus'} ...
);



% %% LSQLIN - node NP
% 
% % Referencia np vektor az optimalizáláshoz az eredeti vektor
% % slack nélküli 14 db csomóponti NP
% x_ref = np;   
% 
% % feltételek együtthatóinak mátrixa
% Aineq = [
%      ntsPTDF_cbco;
%     -ntsPTDF_cbco
% ];
% 
% [~, edgeIndex] = ismember( ...
%     string(CBCO_table.CriticalBranch), string(edges.EdgeName));
% 
% limit = edges.Fmax(edgeIndex) - edges.FRM(edgeIndex);
% 
% % az egyenlőtlenségek jobb oldala
% bineq = [
%     limit;
%     limit
% ];
% 
% % bineq = [
% %      finalram_plus_list;
% %      finalram_minus_list
% % ];
% 
% 
% % egységmátrix, mert min(Cx-xref)^2
% C = eye(numel(x_ref));
% 
% % referencia np a matlabos jelölésnek megfelelően
% d = x_ref;
% 
% % egyenlőségek nincsenek a feltételek között
% % a sum(np)=0 nincs szükség, mert a slack később lesz visszaszámolva
% % matek miatt nem lesz túlterhelődési probléma (PTDF*x szorzatban benne van
% % a slack hatása is)
% % ha np-ben benne van slack, és ptdf-ekben 0-s oszlop,sor, akkor kell
% % sum(np)=0
% Aeq = [];
% beq = [];
% 
% % alsó korlát np-kre -1000
% lb = -1000 * ones(numel(x_ref),1);
% 
% % felső korlát np-kre 1000
% ub =  1000 * ones(numel(x_ref),1);
% 
% options = optimoptions('lsqlin','Display','iter');
% 
% 
% % resnorm: eltérés négyzete
% % primal: x_new (új np)
% % dual: korláthoz tartozik
% % residual:C*x_new - d, vagyis az np változást leíró vektor
% % exitflag: jelzi, hogy sikeres volt-e az optimalizálás.
% [x_new,resnorm,residual,exitflag,output] = lsqlin( ...
%     C,d,Aineq,bineq,Aeq,beq,lb,ub,[],options);
% 
% % slack np beillesztése
% nodeNP_new = zeros(height(nodes),1);
% 
% nodeNP_new(nodes.NodeName ~= slacknodename) = x_new;
% nodeNP_new(slacknodeid) = -sum(x_new);
% 
% nodeNP_result_table = table( ...
%     nodes.NodeName, nodes.NP, nodeNP_new, nodeNP_new - nodes.NP, ...
%     'VariableNames', {'Node','NP_original','NP_new','DeltaNP'})
% 
% 
% % Új np ellenőrzés
% left_plus_new  = ntsPTDF_cbco * x_new;
% left_minus_new = -ntsPTDF_cbco * x_new;
% 
% lsqlin_check_table = table( ...
%     string(CBCO_table.CriticalBranch), ...
%     string(CBCO_table.CriticalOutage), ...
%     left_plus_new, finalram_plus_list, finalram_plus_list - left_plus_new, ...
%     left_minus_new, finalram_minus_list, finalram_minus_list - left_minus_new, ...
%     'VariableNames', {'CB','CO', ...
%     'PTDFxNP_plus','RAM_plus','Difference_plus', ...
%     'PTDFxNP_minus','RAM_minus','Difference_minus'})


%% Validácio es piacszűkítés

result = validation(edges, zonetoslackPTDF, LODF, ...
    zoneNP, zonenames, Ramr);

CBCO_all_table = result.CBCO_table;
reference_table = result.ReferenceTable;
overload_table = result.OverloadedTable;
CBCO_safe_table = result.CBCO_safe_table;

%% ÁBRÁK
edgetable = mergevars(edges,["From","To"],"NewVariableName","EndNodes");
edgetable = movevars(edgetable,"EndNodes","Before",1);

nodetable = renamevars(nodes,"NodeName","Name");
G = digraph(edgetable,nodetable);
% Hálózat alaprajza
G.Nodes.Label = compose("%s\n%+.1f",string(G.Nodes.Name), G.Nodes.NP);
G.Edges.Label = G.Edges.EdgeName;
plot_network(G, slacknodename, 'Mintahálózat');

% Eredeti (referencia) áramlás az éleken
G.Edges.Label = compose("%s: %.1f",G.Edges.EdgeName, G.Edges.Fref);
plot_network(G, slacknodename, 'Referencia áramlás az éleken'); 
% Node-to-node PTDF
mw_ntn = 1;
[~,loc]=ismember(G.Edges.EdgeName,edgetable.EdgeName);
flow_ntn = mw_ntn * nodetonodePTDF(loc,:);
G.Edges.Label = compose("%s: %.1f",G.Edges.EdgeName, flow_ntn);
plot_network(G, slacknodename, sprintf('Áramlás + %g MW node-to-node: %s -> %s', mw_ntn, nodefrom, nodeto)); 

% Node-to-node PTDF után flow
mw_ntn = 1;
flow_ntn = G.Edges.Fref + mw_ntn * nodetonodePTDF(loc,:);
G.Edges.Label = compose("%s: %.1f",G.Edges.EdgeName, flow_ntn);
plot_network(G, slacknodename, sprintf('Áramlás + %g MW node-to-node: %s -> %s', mw_ntn, nodefrom, nodeto));



% Zone-to-slack PTDF egy kiválasztott zónára
[~, loc] = ismember(G.Edges.EdgeName, edges.EdgeName);
G.Edges.Label = compose("%s: %.1f",G.Edges.EdgeName, zonetoslackPTDF(loc,zoneid_plot));
plot_network(G, slacknodename, ... 
    sprintf('Zone-to-slack PTDF: %s -> slack (1 MW)', zoneplot)); 


% Zone-to-slack áramlás egy kiválasztott zónára
mw_zts = 1;
flow_z2s = G.Edges.Fref + mw_zts * zonetoslackPTDF(loc, zoneid_plot);
G.Edges.Label = compose("%s: %.1f",G.Edges.EdgeName, flow_z2s);
plot_network(G, slacknodename, ... 
    sprintf('Áramlás + %g MW zone-to-slack: %s -> slack', ...  
    mw_zts, zoneplot)); 


% Zone-to-zone PTDF adott zónából egy másik, megadott zónába
G.Edges.Label = compose("%s: %.1f",G.Edges.EdgeName, zonetozonePTDF(loc,:));
plot_network(G, slacknodename, sprintf('Zone-to-zone PTDF: %s -> %s (1 MW)', zonefrom, zoneto)) 


% Zone-to-zone áramlás adott zónából egy másik, megadott zónába
mw_ztz = 1;
flow_z2z = G.Edges.Fref + mw_ztz * zonetozonePTDF(loc,:);
G.Edges.Label = compose("%s: %.1f",G.Edges.EdgeName, flow_z2z);
plot_network(G, slacknodename, ... 
    sprintf('Áramlás + %g MW zone-to-zone: %s -> %s', mw_ztz, zonefrom, zoneto)) 

% LODF egy kiválasztott kieső vezetékre
G.Edges.Label = compose("%s: %.1f",G.Edges.EdgeName, LODF(loc,outageid));
plot_network(G, slacknodename, ... 
    sprintf('LODF együtthatók %s kiesése esetén', outageline))

% Kiesés utáni flow
flowaftout = G.Edges.Fref + LODF(loc, outageid) * G.Edges.Fref(G.Edges.EdgeName == outageline);
flowaftout(G.Edges.EdgeName == outageline) = 0;
G.Edges.Label = compose("%s: %.1f",G.Edges.EdgeName, flowaftout);
plot_network(G, slacknodename, ... 
    sprintf('%s kiesése utáni áramlások', outageline))



% Áramlás kereskedelem nélkül
% NP0​=NPref​−GSK⋅NPz​
% F0​=PTDFn2s​⋅NP0
% F0​=Fref​−PTDFz2s​⋅NPz
nodes.NP0 = zeros(height(nodes),1);
nodes.NP0(nodes.NodeName ~= slacknodename) = ...
    nodes.NP(nodes.NodeName ~= slacknodename) - GSK*zoneNP;
nodes.NP0(nodes.NodeName == slacknodename) = -sum(nodes.NP0(nodes.NodeName ~= slacknodename));
G.Nodes.NP0 = nodes.NP0;
G.Nodes.Label = compose("%s\n%+.1f", string(G.Nodes.Name), G.Nodes.NP0);
G.Edges.Label = compose("%s: %.1f",G.Edges.EdgeName, G.Edges.F0);
plot_network(G,slacknodename,'Áramlás kereskedelem nélkül (F0)')


% RAM0+ az egyes vezetékeken F0 esetben
G.Edges.Label = compose("%s: %.1f",G.Edges.EdgeName, G.Edges.RAM0_plus);
plot_network(G,slacknodename,'RAM0+ az egyes vezetékeken')

% RAM0- az egyes vezetékeken F0 esetben
G.Edges.Label = compose("%s: %.1f",G.Edges.EdgeName, G.Edges.RAM0_minus);
plot_network(G,slacknodename,'RAM0- az egyes vezetékeken')

% Virtuális kapacitás pozitív irányban az egyes vezetékeken
G.Nodes.Label = compose("%s\n%+.1f", string(G.Nodes.Name), G.Nodes.NP);
G.Edges.Label = compose("%s: %.1f",G.Edges.EdgeName, G.Edges.AMR_plus);
plot_network(G,slacknodename,'AMR+ az egyes vezetékeken')

% Virtuális kapacitás negatív irányban az egyes vezetékeken
G.Nodes.Label = compose("%s\n%+.1f", string(G.Nodes.Name), G.Nodes.NP);
G.Edges.Label = compose("%s: %.1f",G.Edges.EdgeName, G.Edges.AMR_minus);
plot_network(G,slacknodename,'AMR- az egyes vezetékeken')



function plot_network(G, slacknodename, figTitle)

    figure('Color','w')
    H = plot(G, ...
        'XData', G.Nodes.X, ...
        'YData', G.Nodes.Y, ...
        'EdgeLabel', G.Edges.Label, ...
        "NodeLabel", repelem("",numnodes(G)),...
        'EdgeFontSize',9,...
        'MarkerSize', 40, ...
        'LineWidth', 1.5, ...
        'ArrowSize', 12);        

   
    axis equal
    axis off
    title(figTitle, 'Interpreter','none')

    nodeColors = get_zone_colors(G.Nodes.Zone);

    for node = 1:numnodes(G)
        if G.Nodes.Name(node) == slacknodename
            highlight(H,node,"NodeColor",brighten(nodeColors(node,:),0.5));         
        else
            highlight(H,node,"NodeColor",nodeColors(node,:));            
        end       
    end


    text(G.Nodes.X,G.Nodes.Y,G.Nodes.Label,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',9,'FontWeight','bold');
    
    
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
