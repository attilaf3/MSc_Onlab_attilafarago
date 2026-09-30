function result = validation(edges, zonetoslackPTDF, LODF, zoneNP, zonenames, Ramr)

zoneNP = zoneNP(:);
zonenames = string(zonenames(:));
tol = 1e-6;

%% Teljes CBCO összeállítása
CB = strings(0,1);
CO = strings(0,1);
PTDF_cbco = [];
Fref = [];
F0 = [];
Fmax = [];
FRM = [];

for co = 0:height(edges)
    for cb = 1:height(edges)
        if cb == co
            continue
        end

        if co == 0
            PTDF_row = zonetoslackPTDF(cb,:);
            Fref_row = edges.Fref(cb);
            outage_name = "BaseCase";
        else
            PTDF_row = zonetoslackPTDF(cb,:) + ...
                LODF(cb,co) * zonetoslackPTDF(co,:);
            Fref_row = edges.Fref(cb) + ...
                LODF(cb,co) * edges.Fref(co);
            outage_name = string(edges.EdgeName(co)) + " kiesése";
        end

        CB(end+1,1) = string(edges.EdgeName(cb));
        CO(end+1,1) = outage_name;
        PTDF_cbco(end+1,:) = PTDF_row;
        Fref(end+1,1) = Fref_row;
        F0(end+1,1) = Fref_row - PTDF_row * zoneNP;
        Fmax(end+1,1) = edges.Fmax(cb);
        FRM(end+1,1) = edges.FRM(cb);
    end
end

%% RAM és AMR
minRAM = Ramr * Fmax;
RAM0_plus = Fmax - FRM - F0;
RAM0_minus = Fmax - FRM + F0;
AMR_plus = max(minRAM - RAM0_plus, 0);
AMR_minus = max(minRAM - RAM0_minus, 0);
RAM_plus = RAM0_plus + AMR_plus;
RAM_minus = RAM0_minus + AMR_minus;

CBCO_table = table(CB, CO, Fref, F0, Fmax, FRM, ...
    RAM0_plus, AMR_plus, RAM_plus, RAM0_minus, AMR_minus, RAM_minus, minRAM);

for z = 1:numel(zonenames)
    CBCO_table.(char("PTDF_" + zonenames(z))) = PTDF_cbco(:,z);
end

nCBCO = height(CBCO_table);

%% A referenciaállapot ellenőrzése
refFlow = F0 + PTDF_cbco * zoneNP;
refOverload = max(abs(refFlow) - Fmax, 0);
refFRMExcess = max(abs(refFlow) - (Fmax - FRM), 0);
refMarginPlus = RAM_plus - PTDF_cbco * zoneNP;
refMarginMinus = RAM_minus + PTDF_cbco * zoneNP;

ref_table = table(CB, CO, refFlow, refOverload, ...
    refFRMExcess, refMarginPlus, refMarginMinus);

%% Domain
A_domain = [PTDF_cbco; -PTDF_cbco];
b_domain = [RAM_plus; RAM_minus];

%% LP-s megközelítés
% T csak ábrához és max a vertexhez
T = [eye(3); -ones(1,3)];

% inicializálás: max áramlások
max_plus = zeros(nCBCO,1);
max_minus = zeros(nCBCO,1);
% inicializálás: max áramlást előidéző NP kombináció
worstNP_plus = zeros(nCBCO,4);
worstNP_minus = zeros(nCBCO,4);

options = optimoptions('linprog', 'Display', 'none');

for c = 1:nCBCO
    % -, mert a linprog minimalizál
    % ones(1,4) a sum(NP)
    % NPplus: a maximumot előidéző NP-pont
    % fval_plus: a minimalizált célfüggvény értéke: -PTDFc*NPplus
    [NPplus, fval_plus] = linprog(-PTDF_cbco(c,:).', ...
        A_domain, b_domain, ones(1,4), 0, -inf(4,1), inf(4,1), options);
    [NPminus, fval_minus] = linprog(PTDF_cbco(c,:).', ...
        A_domain, b_domain, ones(1,4), 0, -inf(4,1), inf(4,1), options);

    % teljes max áramlás
    % negatív előjel miatt kivonás
    max_plus(c) = F0(c) - fval_plus;
    max_minus(c) = -F0(c) - fval_minus;
    worstNP_plus(c,:) = NPplus.';
    worstNP_minus(c,:) = NPminus.';
end

%% Vertexes megközelítés
%{
A_reduced = A_domain * T;
polarPoints = unique(bsxfun(@rdivide, A_reduced, b_domain), 'rows');
polarFaces = convhulln(polarPoints);
vertices3 = [];

for j = 1:size(polarFaces,1)
    Qface = polarPoints(polarFaces(j,:),:);
    vertex = Qface \ ones(3,1);
    vertices3(end+1,:) = vertex.';
end

vertices3 = uniquetol(vertices3, 1e-8, 'ByRows', true, 'DataScale', 1);
NP_vertices = vertices3 * T.';
flowsAtVertices = bsxfun(@plus, F0, PTDF_cbco * NP_vertices.');
[max_plus, vertex_plus] = max(flowsAtVertices, [], 2);
[max_minus, vertex_minus] = max(-flowsAtVertices, [], 2);
worstNP_plus = NP_vertices(vertex_plus,:);
worstNP_minus = NP_vertices(vertex_minus,:);

result.VerticesTable = array2table(NP_vertices, ...
    'VariableNames', cellstr("NP_" + zonenames));
result.Vertices = NP_vertices;
result.FlowsAtVertices = flowsAtVertices;
%}

%% Eredmények

% 1-484 pozitív + 484-968 negatív
direction = [repmat("+",nCBCO,1); repmat("-",nCBCO,1)];

% max fellépő áramlás (csak érdekességnek)
maxDirectedFlow = [max_plus; max_minus];

% adott irányban max áramlás (csak érdekességnek)
signedFlow = [max_plus; -max_minus];

% minden sorhoz való NP kombináció hozzárendelése
worstNP = [worstNP_plus; worstNP_minus];

% túlterhelés, [Fmax;Fmax] szükséges, mivel minden CBCO kétszer szerepel
physicalOverload = max(maxDirectedFlow - [Fmax; Fmax], 0);
FRMLevelExcess = max(maxDirectedFlow - [Fmax-FRM; Fmax-FRM], 0);
isOverloaded = FRMLevelExcess > tol;

% CBCO-indexeket kétszeri felsorolása, hogy a pozitív és negatív
% eredménysora ugyanahhoz a CBCO-hoz tartozzon
cbcoRow = [(1:nCBCO).'; (1:nCBCO).'];

results_table = table(cbcoRow, [CB; CB], [CO; CO], direction, ...
    maxDirectedFlow, signedFlow, [Fmax; Fmax], physicalOverload, ...
    FRMLevelExcess, isOverloaded, ...
    'VariableNames', {'CBCORow','CB','CO','Direction','MaxDirectedFlow', ...
    'SignedFlow','Fmax','Overload','FRMLevelExcess','IsOverloaded'});

for z = 1:numel(zonenames)
    results_table.(char("NP_" + zonenames(z))) = worstNP(:,z);
end

result.CBCO_table = CBCO_table;
result.RefTable = ref_table;

% referenciaállaot ellenőrzése
result.RefInDomain = abs(sum(zoneNP)) <= tol && ...
    all([refMarginPlus; refMarginMinus] >= -tol);

result.ResultsTable = results_table;
result.OverloadedTable = sortrows(results_table(isOverloaded,:), 'Overload', 'descend');
result.PTDF_cbco = PTDF_cbco;
result.F0 = F0;
result.Fmax = Fmax;
result.FRM = FRM;
result.A_domain = A_domain;
result.b_domain = b_domain;
result.b_physical = [Fmax-F0; Fmax+F0];
result.T = T;
result.UnsafeDirectionIndices = find(isOverloaded);
result.Tolerance = tol;

%% Piac szűkítése

% biztonságos áramlás
b_physical = [
    Fmax - FRM - F0;
    Fmax - FRM + F0
];
% elemenkénti minimum
b_modified = min(b_domain, b_physical);

result.b_new = b_modified;
result.RAMReduction = b_domain - b_modified;
result.CBCO_safe_table = CBCO_table;
result.CBCO_safe_table.RAM_market_plus = RAM_plus;
result.CBCO_safe_table.RAM_market_minus = RAM_minus;

% oszlop átnevezése
result.CBCO_safe_table.Properties.VariableNames{strcmp(result.CBCO_safe_table.Properties.VariableNames,'AMR_plus')} = 'AMR_market_plus';
result.CBCO_safe_table.Properties.VariableNames{strcmp(result.CBCO_safe_table.Properties.VariableNames,'AMR_minus')} = 'AMR_market_minus';

% új RAM-ok értékeinek beírása
result.CBCO_safe_table.RAM_plus = b_modified(1:nCBCO);
result.CBCO_safe_table.RAM_minus = b_modified(nCBCO+1:end);

% annak ellenőrzése, hogy maradt e még megengedett pont a szűkítés után
% eredményt nem tároljuk
[~, ~, flag] = linprog(zeros(4,1), A_domain, b_modified, ...
    ones(1,4), 0, -inf(4,1), inf(4,1), options);

result.SafeDomainFeasible = flag > 0;
result.SafeDomainVerified = false;
result.PlotAvailable = false;

if flag <= 0
    return
end

%% Az új domain LP-s ellenőrzése
new_max_plus = zeros(nCBCO,1);
new_max_minus = zeros(nCBCO,1);
new_NP_plus = zeros(nCBCO,4);
new_NP_minus = zeros(nCBCO,4);
LP_flags = zeros(nCBCO,2);

for c = 1:nCBCO
    % megkeressük az új pozitív maximumot
    [NPplus, fval_plus, LP_flags(c,1)] = linprog(-PTDF_cbco(c,:).', ...
        A_domain, b_modified, ones(1,4), 0, -inf(4,1), inf(4,1), options);
    % megkeressük az új negatív maximumot
    [NPminus, fval_minus, LP_flags(c,2)] = linprog(PTDF_cbco(c,:).', ...
        A_domain, b_modified, ones(1,4), 0, -inf(4,1), inf(4,1), options);

    if any(LP_flags(c,:) <= 0)
        result.SafeLPFlags = LP_flags;
        return
    end

    new_max_plus(c) = F0(c) - fval_plus;
    new_max_minus(c) = -F0(c) - fval_minus;
    new_NP_plus(c,:) = NPplus.';
    new_NP_minus(c,:) = NPminus.';
end

% egymás alá a két irány eredményeit
new_flow = [new_max_plus; new_max_minus];
new_overload = max(new_flow - [Fmax-FRM; Fmax-FRM], 0);
new_NP = [new_NP_plus; new_NP_minus];

% összehasonlító táblázat
new_table = table(cbcoRow, [CB; CB], [CO; CO], direction, ...
    maxDirectedFlow, new_flow, FRMLevelExcess, new_overload, ...
    'VariableNames', {'CBCORow','CB','CO','Direction', ...
    'OldMaxDirectedFlow','NewMaxDirectedFlow','OldOverload','NewOverload'});

for z = 1:numel(zonenames)
    new_table.(char("NP_" + zonenames(z))) = new_NP(:,z);
end

result.SafeResultsTable = new_table;
result.SafeOverloadedTable = new_table(new_overload > tol,:);
result.SafeDomainVerified = all(new_overload <= tol);
result.SafeLPFlags = LP_flags;

%% Domaincsúcsok az ábrákhoz
vertices_old = domain_vertices(A_domain, b_domain, T, options);
vertices_new = domain_vertices(A_domain, b_modified, T, options);

result.Vertices = vertices_old;
result.SafeVertices = vertices_new;
result.VerticesTable = array2table(vertices_old, ...
    'VariableNames', cellstr("NP_" + zonenames));
result.SafeVerticesTable = array2table(vertices_new, ...
    'VariableNames', cellstr("NP_" + zonenames));

if isempty(vertices_old) || isempty(vertices_new)
    return
end

[faces_old, volume_old] = convhulln(vertices_old(:,1:3));
[faces_new, volume_new] = convhulln(vertices_new(:,1:3));

result.VolumeOld = volume_old;
result.VolumeNew = volume_new;
result.VolumeReductionPercent = 100 * (volume_old-volume_new) / volume_old;

%% 2D ábrák
pairs = nchoosek(1:numel(zonenames),2);
fig = figure('Color','w');

for k = 1:size(pairs,1)
    z1 = pairs(k,1);
    z2 = pairs(k,2);
    points_old = vertices_old(:,[z1 z2]);
    points_new = vertices_new(:,[z1 z2]);
    boundary_old = convhull(points_old(:,1),points_old(:,2));
    boundary_new = convhull(points_new(:,1),points_new(:,2));

    subplot(2,3,k)
    hold on

    fill(points_old(boundary_old,1),points_old(boundary_old,2), ...
        [0.2 0.5 0.9], 'FaceAlpha',0.25, ...
        'EdgeColor',[0.2 0.4 0.7], 'LineWidth',1.5);

    fill(points_new(boundary_new,1),points_new(boundary_new,2), ...
        [0.2 0.7 0.3], 'FaceAlpha',0.45, ...
        'EdgeColor',[0.1 0.4 0.2], 'LineWidth',1.5);

    xlabel(char("NP_" + zonenames(z1)), 'Interpreter','none')
    ylabel(char("NP_" + zonenames(z2)), 'Interpreter','none')
    title(char(zonenames(z1) + " - " + zonenames(z2)))

    axis equal
    grid on
    legend('Kiinduló','Módosított','Location','best')
end

sgtitle('Kiinduló és módosított domain-ek')

result.DomainFigure = fig;
result.PlotAvailable = true;
end

function vertices = domain_vertices(A_domain, b_domain, T, options)

A_reduced = A_domain * T;
row_length = sqrt(sum(A_reduced.^2,2));

[center, ~, flag] = linprog([0;0;0;-1], ...
    [A_reduced row_length], b_domain, [], [], ...
    [-inf;-inf;-inf;0], [], options);

vertices = zeros(0,4);

if flag <= 0 || center(4) <= 1e-8
    return
end

center = center(1:3);
b_shifted = b_domain - A_reduced * center;
polarPoints = unique(bsxfun(@rdivide, A_reduced, b_shifted), 'rows');
polarFaces = convhulln(polarPoints);
vertices3 = [];

for j = 1:size(polarFaces,1)
    Qface = polarPoints(polarFaces(j,:),:);
    vertex = center + (Qface \ ones(3,1));
    vertices3(end+1,:) = vertex.';
end

vertices3 = uniquetol(vertices3, 1e-8, 'ByRows', true, 'DataScale', 1);
vertices = vertices3 * T.';
end