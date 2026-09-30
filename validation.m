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
result.b_physical = [Fmax-FRM-F0; Fmax-FRM+F0];
% result.b_physical = [Fmax-F0; Fmax+F0];
result.T = T;
result.UnsafeDirectionIndices = find(isOverloaded);
result.Tolerance = tol;

%% Piac szűkítése

% Betartandó határ minden CBCO-ban
flow_limit = Fmax - FRM;
% flow_limit = Fmax;

% pozitív:
% F0 + PTDF*NP <= flow_limit
% Ebből: PTDF*NP <= flow_limit - F0
% negatív:
% -PTDF*NP <= flow_limit + F0
b_physical = [
    flow_limit - F0;
    flow_limit + F0
];

% Kezdeti domain-nel való inicializálás
b_modified = b_domain;

% Az első feladatban megtalált maximumok.
% első nCBCO érték pozitív, a második nCBCO negatív irányú
current_max = maxDirectedFlow;

% ebbe, hogy melyik korlátot mennyivel módosítottuk.
% oszlopok: korlát indexe, régi RAM, új RAM, korlátsértés.
history = zeros(0,4);

result.b_physical = b_physical;
result.SafeDomainFeasible = false;
result.SafeDomainVerified = false;
result.PlotAvailable = false;

while true
    % Megnézzük, mennyivel léphető túl a megengedett áramlás.
    % [flow_limit; flow_limit], mert két irányt vizsgálunk, és így
    % egyszerűbb
    excess = current_max - [flow_limit; flow_limit];

    % jelenlegi legnagyobb korlátsértést.
    % selected az A_domain megfelelő sorának indexe
    [largest_excess, selected] = max(excess);

    % Ha már egyik irányban sincs korlátsértés, megállunk, akkor break
    if largest_excess <= tol
        break
    end

    % figyelt korlát jelenlegi jobb oldala.
    old_RAM = b_modified(selected);

    % csak ennek a korlátnak a szigorítása
    new_RAM = min(old_RAM, b_physical(selected));

 
    % numerikus problémák elkerülésére
    % ne menjen a ciklus a végtelenségig
    if old_RAM-new_RAM <= tol
        result.RestrictionFailedIndex = selected;
        result.b_new = b_modified;
        return
    end

    b_modified(selected) = new_RAM;

    % módosítások mentése
    history(end+1,:) = [
        selected, old_RAM, new_RAM, largest_excess
    ];

    % új domainben az új maximum megtalálása
    % mivel a vágás más túlterhelt pontot is kizárhatott
    current_max = zeros(2*nCBCO,1);

    for c = 1:nCBCO
        % pozitív irányú maximum.
        [~, fval_plus, flag_plus] = linprog( ...
            -PTDF_cbco(c,:).', ...
            A_domain, b_modified, ones(1,4), 0, ...
            -inf(4,1), inf(4,1), options);

        % negatív irányú maximum
        [~, fval_minus, flag_minus] = linprog( ...
            PTDF_cbco(c,:).', ...
            A_domain, b_modified, ones(1,4), 0, ...
            -inf(4,1), inf(4,1), options);

        % Ez például üres domain miatt is történhet.
        if flag_plus <= 0 || flag_minus <= 0
            result.RestrictionFailedCBCO = c;
            result.RestrictionLPFlags = [flag_plus flag_minus];
            result.b_new = b_modified;
            return
        end

        % Az F0 állandó áramlás hozzáadása a célfüggvényhez
        current_max(c) = F0(c)-fval_plus;
        current_max(nCBCO+c) = -F0(c)-fval_minus;
    end

    % következő körben már ezekből az új maximumokból kerül ki a következő
    % legnagyobb korlátsértés
end

% végleges RAM-ok és a soronkénti csökkentések.
result.b_new = b_modified;
result.RAMReduction = b_domain-b_modified;

% módosítások sorrendje egy egyszerű táblázatban.
result.RestrictionHistory = array2table(history, ...
    'VariableNames', ...
    {'ConstraintIndex','OldRAM','NewRAM','ExcessBefore'});

% eredeti CBCO-táblából kiindulás
result.CBCO_safe_table = CBCO_table;

% eredeti RAM-ok
result.CBCO_safe_table.RAM_market_plus = RAM_plus;
result.CBCO_safe_table.RAM_market_minus = RAM_minus;

% eredeti AMR-ek
result.CBCO_safe_table.AMR_market_plus = AMR_plus;
result.CBCO_safe_table.AMR_market_minus = AMR_minus;

% a szűkítés utáni RAM
result.CBCO_safe_table.RAM_plus = b_modified(1:nCBCO);
result.CBCO_safe_table.RAM_minus = b_modified(nCBCO+1:end);

% új  AMR = új RAM - eredeti RAM0.
result.CBCO_safe_table.AMR_plus = ...
    b_modified(1:nCBCO)-RAM0_plus;

result.CBCO_safe_table.AMR_minus = ...
    b_modified(nCBCO+1:end)-RAM0_minus;

% ez is csak anomália detektálásához 
% 0 a célfüggvény, mert itt csak megengedett pont keresése
[~, ~, flag] = linprog(zeros(4,1), ...
    A_domain, b_modified, ones(1,4), 0, ...
    -inf(4,1), inf(4,1), options);

result.SafeDomainFeasible = flag > 0;

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
new_overload = max(new_flow - [flow_limit; flow_limit], 0);
% new_overload = max(new_flow - [Fmax Fmax], 0);
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



%% ábrák


pairs = nchoosek(1:4,2);
domain_figures = cell(6,1);
detail_figures = {};

for k = 1:size(pairs,1)
    z1 = pairs(k,1);
    z2 = pairs(k,2);

    remaining = setdiff(1:4,[z1 z2]);
    fixed_zone = remaining(1);
    balance_zone = remaining(2);

    S = zeros(4,2);
    S(z1,1) = 1;
    S(z2,2) = 1;
    S(balance_zone,:) = [-1 -1];

    NP_fixed = zeros(4,1);
    NP_fixed(fixed_zone) = zoneNP(fixed_zone);
    NP_fixed(balance_zone) = -zoneNP(fixed_zone);

    % Az összes korlát átírása két változóra.
    A2 = A_domain*S;
    offset = A_domain*NP_fixed;

    b2_old = b_domain-offset;
    b2_new = b_modified-offset;
    b2 = [b2_old b2_new];

    points = cell(2,1);

    % A régi és az új metszet vertexei.
    for domain = 1:2
        points{domain} = zeros(0,2);

        row_length = sqrt(sum(A2.^2,2));

        [center, ~, flag] = linprog([0;0;-1], ...
            [A2 row_length], b2(:,domain), [], [], ...
            [-inf;-inf;0], [], options);

        if flag <= 0
            continue
        end

        if center(3) <= tol
            continue
        end

        center = center(1:2);
        b_shifted = b2(:,domain)-A2*center;

        polarPoints = unique( ...
            bsxfun(@rdivide,A2,b_shifted),'rows');

        polarBoundary = convhull( ...
            polarPoints(:,1),polarPoints(:,2));

        vertices2 = [];

        for j = 1:numel(polarBoundary)-1
            Q = polarPoints(polarBoundary(j:j+1),:);
            vertex = center+(Q \ ones(2,1));
            vertices2(end+1,:) = vertex.';
        end

        vertices2 = uniquetol(vertices2,1e-8, ...
            'ByRows',true,'DataScale',1);

        boundary = convhull(vertices2(:,1),vertices2(:,2));
        points{domain} = vertices2(boundary(1:end-1),:);
    end

    p_old = points{1};
    p_new = points{2};

    if isempty(p_old)
        continue
    end

    low = min(p_old,[],1);
    high = max(p_old,[],1);
    padding = 0.25*(high-low);

    limits = [
        low(1)-padding(1), high(1)+padding(1), ...
        low(2)-padding(2), high(2)+padding(2)
    ];

    % A régi vagy új határt alkotó szigorított korlátok.
    cuts = [];

    for c = 1:size(A2,1)
        if b2_old(c)-b2_new(c) <= tol
            continue
        end

        old_on_line = abs(p_old*A2(c,:).'-b2_old(c)) < 10*tol;
        new_on_line = false(size(p_new,1),1);

        if ~isempty(p_new)
            new_on_line = ...
                abs(p_new*A2(c,:).'-b2_new(c)) < 10*tol;
        end

        if nnz(old_on_line) >= 2 || nnz(new_on_line) >= 2
            cuts(end+1,1) = c;
        end
    end

    % A 0 az áttekintő ábra.
    % Az A–B metszetnél külön részletes ábrák is készülnek.
    plot_rows = 0;

    if k == 1
        plot_rows = [0; cuts];
    end

    for picture = 1:numel(plot_rows)
        selected = plot_rows(picture);

        fig = figure('Color','w');
        hold on

        if selected == 0
            domain_figures{k} = fig;
            rows_to_draw = cuts;
        else
            detail_figures{end+1} = fig;
            rows_to_draw = selected;
        end

        % Eredeti domain, pirosas háttérrel.
        fill(p_old(:,1),p_old(:,2),[1 0.8 0.75], ...
            'EdgeColor',[0.2 0.4 0.7], ...
            'LineWidth',2,'HandleVisibility','off');

        % Az összes CBCO-val szűkített domain.
        if ~isempty(p_new)
            fill(p_new(:,1),p_new(:,2),[0.75 0.9 0.8], ...
                'EdgeColor',[0.1 0.4 0.2], ...
                'LineWidth',2,'HandleVisibility','off');
        end

        for c = rows_to_draw.'
            PTDF_row = A2(c,:);
            limit = b2_new(c);

            if PTDF_row*PTDF_row.' < 1e-12
                continue
            end

            % Régi és eltolt határoló egyenes.
            if abs(PTDF_row(2)) > 1e-12
                x = linspace(limits(1),limits(2),100);

                y_old = (b2_old(c)-PTDF_row(1)*x) ...
                    / PTDF_row(2);

                y_new = (limit-PTDF_row(1)*x) ...
                    / PTDF_row(2);

                plot(x,y_old,'b--','LineWidth',1.3, ...
                    'HandleVisibility','off');

                plot(x,y_new,'r-','LineWidth',1.5, ...
                    'HandleVisibility','off');
            else
                x_old = b2_old(c)/PTDF_row(1);
                x_new = limit/PTDF_row(1);

                plot([x_old x_old],limits(3:4),'b--', ...
                    'LineWidth',1.3,'HandleVisibility','off');

                plot([x_new x_new],limits(3:4),'r-', ...
                    'LineWidth',1.5,'HandleVisibility','off');
            end

            number = 0;

            % A régi élek metszése az új korláttal.
            for j = 1:size(p_old,1)
                v1 = p_old(j,:);
                v2 = p_old(mod(j,size(p_old,1))+1,:);

                d1 = PTDF_row*v1.'-limit;
                d2 = PTDF_row*v2.'-limit;

                if ~((d1 < -tol && d2 > tol) || ...
                        (d2 < -tol && d1 > tol))
                    continue
                end

                % V1 biztonságos, V2 sértő az adott korlátra.
                if d1 > 0
                    temporary = v1;
                    v1 = v2;
                    v2 = temporary;
                end

                t = (limit-PTDF_row*v1.') / ...
                    (PTDF_row*(v2-v1).');

                w = v1+t*(v2-v1);
                number = number+1;

                % Merőleges eltolási nyíl a régi határtól.
                start = w + (b2_old(c)-limit) / ...
                    (PTDF_row*PTDF_row.') * PTDF_row;

                movement = w-start;

                quiver(start(1),start(2), ...
                    movement(1),movement(2),0, ...
                    'Color','c','LineWidth',1.5, ...
                    'MaxHeadSize',0.4, ...
                    'HandleVisibility','off');

                % A részletes ábrán az él és a pontok jelölése.
                if selected ~= 0
                    plot([v1(1) v2(1)],[v1(2) v2(2)], ...
                        'Color',[0.9 0.5 0.1], ...
                        'LineWidth',3,'HandleVisibility','off');

                    plot(v1(1),v1(2),'bo', ...
                        'MarkerFaceColor','b','MarkerSize',7, ...
                        'HandleVisibility','off');

                    plot(v2(1),v2(2),'rx', ...
                        'LineWidth',2,'MarkerSize',10, ...
                        'HandleVisibility','off');

                    plot(w(1),w(2),'mo', ...
                        'MarkerFaceColor','m','MarkerSize',8, ...
                        'HandleVisibility','off');

                    text(v1(1),v1(2), ...
                        sprintf('  V1(%d)',number),'Color','b');

                    text(v2(1),v2(2), ...
                        sprintf('  V2(%d): sértő',number), ...
                        'Color','r');

                    text(w(1),w(2), ...
                        sprintf('  W%d; t = %.3f',number,t), ...
                        'Color','m');
                end
            end
        end

        % Az áttekintő ábrákon a régi és új vertexek.
        if selected == 0
            old_is_safe = all( ...
                A2*p_old.' <= b2_new+tol,1).';

            plot(p_old(old_is_safe,1), ...
                p_old(old_is_safe,2),'bo', ...
                'MarkerFaceColor','b', ...
                'HandleVisibility','off');

            plot(p_old(~old_is_safe,1), ...
                p_old(~old_is_safe,2),'rx', ...
                'LineWidth',1.5,'MarkerSize',8, ...
                'HandleVisibility','off');

            if ~isempty(p_new)
                plot(p_new(:,1),p_new(:,2),'mo', ...
                    'MarkerFaceColor','m', ...
                    'HandleVisibility','off');
            end
        end

        h1 = plot(nan,nan,'b--');
        h2 = plot(nan,nan,'r-');
        h3 = plot(nan,nan,'s', ...
            'Color',[0.1 0.4 0.2], ...
            'MarkerFaceColor',[0.75 0.9 0.8]);

        legend([h1 h2 h3], ...
            {'Eredeti korlát','Eltolt korlát', ...
            'Összes CBCO-val szűkített domain'}, ...
            'Location','best');

        xlabel(char("NP_"+zonenames(z1)),'Interpreter','none')
        ylabel(char("NP_"+zonenames(z2)),'Interpreter','none')

        if selected == 0
            title(sprintf('%s–%s; NP_%s = %.4f', ...
                zonenames(z1),zonenames(z2), ...
                zonenames(fixed_zone),zoneNP(fixed_zone)), ...
                'Interpreter','none');
        else
            if selected <= nCBCO
                cbco = selected;
                direction_text = '+';
            else
                cbco = selected-nCBCO;
                direction_text = '-';
            end

            title(sprintf('%s / %s; %s irány; NP_%s = %.4f', ...
                CB(cbco),CO(cbco),direction_text, ...
                zonenames(fixed_zone),zoneNP(fixed_zone)), ...
                'Interpreter','none');
        end

        axis equal
        axis(limits)
        grid on
    end
end

result.DomainFigures = domain_figures;
result.DetailFigures = detail_figures;
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