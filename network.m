function [nodes, edges, slacknodename, slacknodeid, noslacknodesname] = build_network()

% Node table
nodesname = ["N1"; "N2"; "N3"; "N4"; "N5"; "N6"; "N7"; "N8"; "N9"; "N10"; "N11"; "N12"; "N13"; "N14"; "N15"];
nodesnp = [190; 70; -80; -60; 200; -40; -130; 0; 120; 130; -150; 60; -170; 40; -120];
nodesnp(end) = -(sum(nodesnp) - nodesnp(end));
nodesx  = 2 * [0.0; 0.0; 0.0; 0.0; 2.2; 3.0; 3.0; 4.2; 4.7; 6.8; 7.8; 8.2; 5.8; 6.0; 8.0];
nodesy  = 2 * [5.0; 3.7; 2.2; 0.4; 3.4; 2.0; 0.6; 5.1; 3.5; 4.0; 4.7; 2.8; 2.2; 1.0; 0.7];
nodezone = ["A";"A";"B";"B";"A";"B";"B";"C";"C";"C";"C";"D";"D";"D";"D"];

nodes = table(nodesname, nodesnp, nodesx, nodesy, nodezone,'VariableNames', {'NodeName', 'NP', 'X', 'Y', 'Zone'});

% Edge table
edgesname = ["L1"; "L2"; "L3"; "L4"; "L5"; "L6"; "L7"; "L8"; "L9"; "L10"; "L11"; "L12"; "L13"; "L14"; "L15"; "L16"; "L17"; "L18"; "L19"; "L20"; "L21"; "L22"];
edgeslength = [150; 125; 175; 160; 100; 170; 140; 115; 145; 210; 180; 215; 135; 130; 230; 155; 150; 160; 165; 155; 145; 235];
edgesohm = 0.4 * edgeslength;
edgesfrom = ["N1"; "N2"; "N3"; "N4"; "N6"; "N3"; "N3"; "N2"; "N5"; "N1"; "N8"; "N8"; "N9"; "N10"; "N11"; "N10"; "N9"; "N13"; "N6"; "N7"; "N14"; "N12"];
edgesto   = ["N2"; "N3"; "N4"; "N7"; "N7"; "N6"; "N5"; "N5"; "N9"; "N8"; "N9"; "N11"; "N10"; "N11"; "N12"; "N12"; "N13"; "N12"; "N13"; "N14"; "N15"; "N15"];
edgesfmax = 170 * ones(numel(edgesname),1);  
edgesfrm  = 5 * ones(numel(edgesname),1);    

edges = table(edgesname, edgesohm, edgesfrom, edgesto, edgesfmax, edgesfrm, 'VariableNames', {'EdgeName', 'Ohm', 'From', 'To', 'Fmax', 'FRM'});

% Slack csomópont
slacknodename = "N10";
slacknodeid = find(nodes.NodeName == slacknodename);
sumnp = sum(nodes.NP) - nodes{slacknodeid, "NP"};
slacknodenp = -sumnp;
nodes.NP(slacknodeid) = slacknodenp;

noslacknodesname = nodes.NodeName;
noslacknodesname(slacknodeid) = [];

end

%[appendix]
%---
