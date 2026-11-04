%[text] # Teljes validáció IVA\-val
%[text] A függvény minden CBCO\-ra megkeresi a piaci domainben lehetséges legnagyobb pozitív és negatív irányú áramlást. Ha FRM\-határsértés lehetséges, a legnagyobb sértéshez tartozó korlát RAM\-ját csökkenti, majd újraszámítja a maximumokat.
%[text] A PTDF és a F0 a validáció alatt változatlan. A piacszűkítés a RAM\-ok csökkentésével történik.
%[text]{"align":"center"} $F=F\_0\+PTDF\_{CBCO}NP$
%[text]{"align":"center"} $RAM\_{validated}^{\\pm}=RAM0^{\\pm}\+AMR\_{market}^{\\pm}\-IVA^{\\pm}$
%[text] A zoneNP(:) és a zonenames(:) oszlopvektorrá alakítja a bemeneteket. A tol csak a numerikus problémák elkerüléséhez.
%[text] ## Fontosabb változók
%[text:table]
%[text] | Változó | Jelentés |
%[text] | --- | --- |
%[text] | `CB`, `CO` | Vizsgált vezeték és kiesés |
%[text] | `zoneNP` | A refprog zónás nettó pozíciója |
%[text] | `PTDF_cbco` | Zónás érzékenységek minden CBCO\-ra |
%[text] | `Fref`, `F0` | Referenciaáramlás és zero\-balance áramlás |
%[text] | `Fmax`, `FRM` | Fizikai határ és megbízhatósági tartalék |
%[text] | `RAM0_plus`, `RAM0_minus` | AMR nélküli kiinduló RAM |
%[text] | `AMR_plus`, `AMR_minus` | Helyi változóként az eredeti AMR |
%[text] | `A_domain`, `b_domain` | Az eredeti piaci domain |
%[text] | `b_modified` | A szűkített domain jobb oldala |
%[text] | `maxDirectedFlow` | Az eredeti domain irányonkénti maximumai |
%[text] | `current_max` | Az aktuális domain irányonkénti maximumai |
%[text] | `selected` | A szigorítandó irány indexe |
%[text] | `IVA` | RAM\-csökkentés az eredeti domainhez képest |
%[text] | `IVA_applied` | Történt\-e közvetlen korrekció az adott irányon |
%[text] | `result.CorrectedTable` | A közvetlenül módosított irányok táblája |
%[text] | `result.RestrictionHistory` | A módosítások sorrendje |
%[text:table]
%[text] A mintahálózatban 22 alapeseti és 22 × 21 kieséses sor készül, összesen 484 CBCO. A két áramlási irány miatt 968 irányított korlátot vizsgálunk.
function result = validation(edges, zonetoslackPTDF, LODF, zoneNP, zonenames, Ramr)

zoneNP = zoneNP(:);
zonenames = string(zonenames(:));
tol = 1e-6;

%%\
%[text] ## Teljes CBCO összeállítása
%[text] A külső ciklus a kontingenciákon, a belső ciklus a vizsgált vezetékeken halad végig. A co = 0 az alapeset, a többi co érték az adott vezeték kiesését jelenti.
%[text] A cb == co esetet kimarad evidensen.
%[text] Kieséskor a referenciaáramlást és a PTDF\-et a LODF segítségével módosítjuk.
%[text]{"align":"center"} $PTDF\_{cb,co}=PTDF\_{cb}\+LODF(cb,co)PTDF\_{co}$
%[text]{"align":"center"} $F\_{ref,cb,co}=F\_{ref,cb}\+LODF(cb,co)F\_{ref,co}$
%[text] A F0\-t a referenciaállapotból számítjuk vissza.
%[text]{"align":"center"} $F\_{0,cb,co}=F\_{ref,cb,co}\-PTDF\_{cb,co}NP\_{ref}$
%[text] Az end\+1 új sort fűz a táblázatba. Ugyanez a sorindex minden vektorban és a PTDF\-mátrixban ugyanazt a CBCO\-t azonosítja.
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

%%\
%[text] ## RAM és AMR
%[text] Az AMR nélküli RAM\-ok:
%[text]{"align":"center"} $RAM0^\+=F\_{max}\-FRM\-F\_0$
%[text]{"align":"center"} $RAM0^\-=F\_{max}\-FRM\+F\_0$
%[text] A 70%\-os minimum\-RAM követelménnyel kiegészítve kapott AMR:
%[text]{"align":"center"} $minRAM=Ramr\\,F\_{max}$
%[text]{"align":"center"} $AMR^{\\pm}=\\max(minRAM\-RAM0^{\\pm},0)$
%[text] A validáció előtti piaci RAM:
%[text]{"align":"center"} $RAM\_{market}^{\\pm}=RAM0^{\\pm}\+AMR^{\\pm}$
%[text] A CBCO\_table az azonos CBCO\-hoz tartozó adatokat rendezi egy sorba. A következő ciklus zónánként hozzáadja a PTDF\-oszlopokat.
%[text] A nCBCO a táblázat sorainak száma, amelyet a későbbi ciklusok és vektorok méretezéséhez használunk.
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

%%\
%[text] ## Referenciaállapot ellenőrzése
%[text] A refFlow a refprog nettó pozícióját helyettesíti az áramlási összefüggésbe.
%[text]{"align":"center"} $refFlow=F\_0\+PTDF\_{CBCO}NP\_{ref}=F\_{ref}$
%[text] A refOverload a fizikai Fmax túllépését mutatja. A refFRMExcess az FRM\-mel csökkentett határ sértését jelzi.
%[text]{"align":"center"} $refOverload=\\max(|refFlow|\-F\_{max},0)$
%[text]{"align":"center"} $refFRMExcess=\\max(|refFlow|\-(F\_{max}\-FRM),0)$
%[text] A refMarginPlus és refMarginMinus a refprog piaci tartalékai. Ezek a megfelelő jobb oldal és bal oldal különbségei.
%[text]{"align":"center"} $refMarginPlus=RAM^\+\-PTDF\_{CBCO}NP\_{ref}$
%[text]{"align":"center"} $refMarginMinus=RAM^\-\+PTDF\_{CBCO}NP\_{ref}$
%[text] Nemnegatív tartalék esetén a megfelelő piaci feltétel teljesül. A ref\_table ezeket az eredményeket tárolja.
refFlow = F0 + PTDF_cbco * zoneNP;
refOverload = max(abs(refFlow) - Fmax, 0);
refFRMExcess = max(abs(refFlow) - (Fmax - FRM), 0);
refMarginPlus = RAM_plus - PTDF_cbco * zoneNP;
refMarginMinus = RAM_minus + PTDF_cbco * zoneNP;

ref_table = table(CB, CO, refFlow, refOverload, ...
    refFRMExcess, refMarginPlus, refMarginMinus);

%%\
%[text] ## Domain
%[text] A pozitív és negatív irányú feltételeket közös kisebb\-egyenlő alakban:
%[text]{"align":"center"} $PTDF\_{CBCO}NP\\leq RAM^\+$
%[text]{"align":"center"} $\-PTDF\_{CBCO}NP\\leq RAM^\-$
%[text] A mátrixok egymás alá helyezésével kapjuk az A\_domain és b\_domain változókat.
%[text]{"align":"center"} $A\_{domain}=\\left[\\matrix{PTDF\_{CBCO}\\cr \-PTDF\_{CBCO}}\\right]$
%[text]{"align":"center"} $b\_{domain}=\\left[\\matrix{RAM\_{plus}\\cr RAM\_{minus}}\\right]$
%[text] Az első nCBCO sor pozitív, a második nCBCO sor negatív irányú. Ezt a sorrendet a maximumoknál és az IVA\-nál is megtartjuk.
%[text] Az NP\-k összegét minden LP\-ben külön egyenlőségként nullára állítjuk.
%[text]{"align":"center"} $A\_{domain}NP\\leq b\_{domain},\\qquad \\sum\_z NP\_z=0$
A_domain = [PTDF_cbco; -PTDF_cbco];
b_domain = [RAM_plus; RAM_minus];

%%\
%[text] ## LP\-s megközelítés
%[text] Az optimalizálási változó a négy zóna nettó pozíciója. Minden CBCO\-ra külön keressük a pozitív és negatív irányú legnagyobb áramlást a domainben.
%[text]{"align":"center"} $max\_{plus}^{(c)}=F\_0^{(c)}\+\\max\_{NP}(PTDF\_{CBCO}(c,:)NP)$
%[text]{"align":"center"} $max\_{minus}^{(c)}=\-F\_0^{(c)}\+\\max\_{NP}(\-PTDF\_{CBCO}(c,:)NP)$
%[text] A linprog minimalizál, ezért a pozitív maximumhoz a PTDF\-sor negatívját adjuk meg célfüggvényként. A transzponálás oszlopvektort készít az együtthatókból.
%[text] A F0 az adott CBCO\-ra állandó, ezért az LP után adjuk hozzá az eredményhez.
%[text]{"align":"center"} $max\_{plus}^{(c)}=F\_0^{(c)}\-fval\_{plus}$
%[text]{"align":"center"} $max\_{minus}^{(c)}=\-F\_0^{(c)}\-fval\_{minus}$
%[text] Az ones(1,4) és a nulla jobb oldal biztosítja a zónás egyensúlyt. A \-inf és inf azt jelenti, hogy külön zónánkénti alsó és felső korlátot itt nem adunk meg.
%[text] A worstNP\_plus és worstNP\_minus a maximumokat előidéző piaci pontokat tárolja. A transzponálás után ezeket soronként írjuk a mátrixokba. A Display = none kikapcsolja a solver iterációs kiírását.
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

%%\
%[text] ## Eredmények
%[text] A direction a sorok irányát jelöli. A repmat a pozitív és negatív jelölést nCBCO alkalommal ismétli meg.
%[text] A maxDirectedFlow az eredeti domain irányonkénti maximumait tartalmazza, és később a piacszűkítés kezdőadata is lesz.
%[text] A negatív irányú rész a mínusz F legnagyobb értékét tárolja.
%[text] A FRMLevelExcess az FRM\-mel csökkentett áramlási határ túllépése.
%[text]{"align":"center"} $FRMLevelExcess=\\max(maxDirectedFlow\-(F\_{max}\-FRM),0)$
%[text] A határértékeket kétszer írjuk egymás alá, mert minden CBCO két irányban szerepel.
%[text] Az isOverloaded a toleranciánál nagyobb sértéseket jelöli. A logikai indexelés ezekből választja ki az OverloadedTable sorait, a sortrows pedig csökkenő határsértés szerint rendezi őket.
%[text] A cbcoRow összekapcsolja a két irányított sort az eredeti CBCO\-sorral. Az NP\-oszlopok a maximumokhoz tartozó piaci pontokat mutatják.
%[text] A RefInDomain akkor igaz, ha a refprog kiegyensúlyozott, és minden eredeti piaci feltételt teljesít.
%[text] A result struktúrában tárolt adatok a függvényen kívül is elérhetők, például result.ResultsTable néven.
% 1-484 pozitív + 484-968 negatív
direction = [repmat("+",nCBCO,1); repmat("-",nCBCO,1)];

% max fellépő áramlás (csak érdekességnek)
maxDirectedFlow = [max_plus; max_minus];

% minden sorhoz való NP kombináció hozzárendelése
worstNP = [worstNP_plus; worstNP_minus];

% túlterhelés, [Fmax;Fmax] szükséges, mivel minden CBCO kétszer szerepel
FRMLevelExcess = max(maxDirectedFlow - [Fmax-FRM; Fmax-FRM], 0);
isOverloaded = FRMLevelExcess > tol;

% CBCO-indexeket kétszeri felsorolása, hogy a pozitív és negatív
% eredménysora ugyanahhoz a CBCO-hoz tartozzon
cbcoRow = [(1:nCBCO).'; (1:nCBCO).'];

results_table = table(cbcoRow, [CB; CB], [CO; CO], direction, ...
    maxDirectedFlow, [Fmax; Fmax], ...
    FRMLevelExcess, isOverloaded, ...
    'VariableNames', {'CBCORow','CB','CO','Direction','MaxDirectedFlow', ...
    'Fmax','FRMLevelExcess','IsOverloaded'});

for z = 1:numel(zonenames)
    results_table.(char("NP_" + zonenames(z))) = worstNP(:,z);
end

result.CBCO_table = CBCO_table;
result.RefTable = ref_table;

% referenciaállaot ellenőrzése
result.RefInDomain = abs(sum(zoneNP)) <= tol && ...
    all([refMarginPlus; refMarginMinus] >= -tol);

result.ResultsTable = results_table;
result.OverloadedTable = sortrows(results_table(isOverloaded,:), 'FRMLevelExcess', 'descend');
result.PTDF_cbco = PTDF_cbco;
result.F0 = F0;
result.Fmax = Fmax;
result.FRM = FRM;
result.A_domain = A_domain;
result.b_domain = b_domain;
result.b_physical = [Fmax-FRM-F0; Fmax-FRM+F0];
result.UnsafeDirectionIndices = find(isOverloaded);
result.Tolerance = tol;

%%\
%[text] ## Piac szűkítése
%[text] A piacszűkítés célja az FRM\-mel csökkentett áramlási határ betartása a teljes megengedett domainben.
%[text]{"align":"center"} $flow\_{limit}=F\_{max}\-FRM$
%[text] Az áramlási feltételekből kapott fizikai RAM\-vektor:
%[text]{"align":"center"} $b\_{physical}=\\left[\\matrix{flow\_{limit}\-F\_0\\cr flow\_{limit}\+F\_0}\\right]$
%[text] A b\_modified kezdetben az eredeti RAM\-vektor másolata. Ebben vezetjük a szigorításokat, miközben a b\_domain megőrzi az eredeti értékeket.
%[text] Az IVA kezdetben nulla. Az eredeti piaci RAM\-hoz képesti teljes csökkentést tárolja.
%[text]{"align":"center"} $IVA=b\_{domain}\-b\_{modified}$
%[text] A current\_max kezdeti értékeit már ismerjük az eredeti domain LP\-s vizsgálatából.
%[text] A history négy oszlopban tárolja a módosítások adatait: korlátindex, régi RAM, új RAM és módosítás előtti határsértés.
% Betartandó határ minden CBCO-ban
flow_limit = Fmax - FRM;

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

% IVA: az eredeti piaci RAM-hoz képesti csökkentés.
% Kezdetben minden irányon nulla.
IVA = zeros(size(b_domain));

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


%[text] ## A következő korrekció kiválasztása
%[text] Az excess az aktuális maximum és a megengedett áramlás különbsége. A max megadja a legnagyobb sértést és annak selected indexét.
%[text]{"align":"center"} $excess=current\_{max}\-\\left[\\matrix{flow\_{limit}\\cr flow\_{limit}}\\right]$
%[text] A selected az A\_domain és b\_modified ugyanazon irányított korlátjára mutat. Ha a legnagyobb sértés sem haladja meg a toleranciát, a break befejezi a ciklust.
%[text] A kiválasztott RAM\-ot a fizikai értékre csökkentjük. A min biztosítja, hogy a művelet ne növelje a jelenlegi RAM\-ot.
%[text]{"align":"center"} $new\_{RAM}=\\min(old\_{RAM},b\_{physical,selected})$
%[text] A hozzá tartozó IVA a tényleges RAM\-csökkentés.
%[text]{"align":"center"} $IVA\_{selected}=b\_{domain,selected}\-new\_{RAM}$
%[text] Ebben az algoritmusban a kiválasztott irány a RAM0\-ra kerül vissza, így azon az irányon a teljes eredeti AMR\-kiegészítést visszavonjuk.
%[text] Minden korrekció után újraszámítjuk az összes maximumot a b\_modified használatával. Egy korlát szigorítása más CBCO\-k veszélyes piaci pontjait is kizárhatja.
%[text] Az old\_RAM\-new\_RAM feltétel az előrelépés hiányát jelzi. A solverállapothoz tartozó return pedig megszakítja a függvényt, ha valamelyik LP nem ad megfelelő megoldást.
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

    % Az adott irány teljes IVA-ja az eredeti piaci RAM-hoz képest.
    IVA(selected) = b_domain(selected)-b_modified(selected);

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

% végleges RAM-ok.
result.b_new = b_modified;

% módosítások sorrendje egy egyszerű táblázatban.
result.RestrictionHistory = array2table(history, ...
    'VariableNames', ...
    {'ConstraintIndex','OldRAM','NewRAM','ExcessBefore'});

%%\
%[text] ## IVA\-értékek és korrekciójelzők
%[text] Az IVA\-vektor első nCBCO eleme pozitív, a többi negatív irányú.
%[text] Az IVA\_applied logikai vektor a toleranciánál nagyobb közvetlen RAM\-csökkentéseket jelöli.
%[text]{"align":"center"} $IVA\_{applied,i}=(IVA\_i>tol)$
%[text] A find az igaz értékek helyét adja vissza. Az IVAIndices irányított korlátindexeket, a CorrectedCBCOIndices pedig CBCO\-indexeket tartalmaz.
%[text] A logikai vagy művelet miatt egy CBCO akkor is módosítottnak számít, ha csak az egyik irányán történt IVA.
%[text] Az IVATable minden irány eredeti RAM\-ját, IVA\-ját és validált RAM\-ját megőrzi. A CorrectedTable logikai indexeléssel csak a módosított sorokat választja ki.
%[text] Az InitiallyUnsafe a kezdeti határsértést, az IVA\_applied a közvetlen korrekciót jelzi. Egy kezdetben sértő irány saját IVA nélkül is biztonságossá válhat más korlátok szigorítása miatt.
IVA_plus = IVA(1:nCBCO);
IVA_minus = IVA(nCBCO+1:end);

IVA_applied = IVA > tol;
IVA_applied_plus = IVA_plus > tol;
IVA_applied_minus = IVA_minus > tol;

result.IVA = IVA;
result.IVA_plus = IVA_plus;
result.IVA_minus = IVA_minus;

result.IVA_applied = IVA_applied;
result.IVA_applied_plus = IVA_applied_plus;
result.IVA_applied_minus = IVA_applied_minus;

% Irányított korlátindexek.
result.IVAIndices = find(IVA_applied);

% CBCO-indexek: legalább az egyik irányon történt IVA.
result.CorrectedCBCOIndices = find( ...
    IVA_applied_plus | IVA_applied_minus);

% Minden irány eredeti RAM-ja, IVA-ja és új RAM-ja.
IVA_table = table( ...
    (1:2*nCBCO).', cbcoRow, [CB;CB], [CO;CO], direction, ...
    b_domain, IVA, b_modified, IVA_applied, isOverloaded, ...
    'VariableNames', {'ConstraintIndex','CBCORow','CB','CO', ...
    'Direction','RAM_market','IVA','RAM_validated', ...
    'IVA_applied','InitiallyUnsafe'});

result.IVATable = IVA_table;

% Csak a közvetlenül IVA-val módosított irányok.
result.CorrectedTable = IVA_table(IVA_applied,:);


%[text] ## A korrekciók sorrendje
%[text] A RestrictionHistory eredeti korlátindexeihez hozzárendeljük a vizsgált vezetéket, a kiesést és az irányt. Ehhez a pozitív és negatív irányok sorrendjében megismételjük a CB és CO neveket.
%[text] Az Iteration a módosítások sorrendjét mutatja. Az IVA\_step az adott lépés RAM\-csökkentése.
%[text]{"align":"center"} $IVA\_{step}=OldRAM\-NewRAM$
%[text] Az ExcessBefore a módosítás előtti áramlási határsértés, az IVA\_total pedig az adott irány eredeti RAM\-jához képesti végső csökkentése.
% A meglévő history táblázat kiegészítése.
all_CB = [CB;CB];
all_CO = [CO;CO];

selected_rows = result.RestrictionHistory.ConstraintIndex;

result.RestrictionHistory.Iteration = ...
    (1:height(result.RestrictionHistory)).';

result.RestrictionHistory.CB = all_CB(selected_rows);
result.RestrictionHistory.CO = all_CO(selected_rows);
result.RestrictionHistory.Direction = direction(selected_rows);

result.RestrictionHistory.IVA_step = ...
    result.RestrictionHistory.OldRAM ...
    -result.RestrictionHistory.NewRAM;

result.RestrictionHistory.IVA_total = IVA(selected_rows);

%%\
%[text] ## Validált CBCO\-tábla
%[text] A CBCO\_safe\_table az eredeti CBCO\-tábla másolatából indul. A RAM\_market és AMR\_market oszlopok megőrzik a validáció előtti értékeket.
%[text] A RAM\_plus és RAM\_minus oszlopokba az új RAM kerül.
%[text]{"align":"center"} $RAM\_{validated}^{\\pm}=RAM0^{\\pm}\+AMR\_{market}^{\\pm}\-IVA^{\\pm}$
%[text] A meglévő kód az AMR\_plus és AMR\_minus oszlopokban a RAM0 felett megmaradó kiegészítést tárolja.
%[text]{"align":"center"} $AMR\_{remaining}^{\\pm}=RAM\_{validated}^{\\pm}\-RAM0^{\\pm}=AMR\_{market}^{\\pm}\-IVA^{\\pm}$
%[text] Az IVA\-értékeket és a korrekciójelzőket külön oszlopokhoz adjuk. A CBCO\-szintű IVA\_applied akkor igaz, ha legalább az egyik irány módosult.
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

% IVA és irányonkénti korrekciójelzők.
result.CBCO_safe_table.IVA_plus = IVA_plus;
result.CBCO_safe_table.IVA_minus = IVA_minus;

result.CBCO_safe_table.IVA_applied_plus = IVA_applied_plus;
result.CBCO_safe_table.IVA_applied_minus = IVA_applied_minus;

% CBCO-szinten igaz, ha legalább az egyik irány módosult.
result.CBCO_safe_table.IVA_applied = ...
    IVA_applied_plus | IVA_applied_minus;

% Az eredeti eredménytáblában is látszódjon a korrekció.
result.ResultsTable.IVA = IVA;
result.ResultsTable.IVA_applied = IVA_applied;


%[text] ## Megengedett pont keresése
%[text] A nulla célfüggvényű linprog csak egy megengedett piaci pontot keres. Itt minden megengedett pont célfüggvényértéke ugyanaz.
%[text] Ezzel azt vizsgáljuk, hogy a szűkítés után maradt\-e megengedett pont a domainben.
%[text]{"align":"center"} $A\_{domain}NP\\leq b\_{modified},\\qquad \\sum\_z NP\_z=0$
%[text] A SafeDomainFeasible a keresés eredményét jelzi. A nem üres domain még nem önmagában igazolja minden áramlási korlát teljesülését; ezt a következő maximumkeresés vizsgálja.
% csak anomália detektálásához 
% 0 a célfüggvény, mert itt csak megengedett pont keresése
[~, ~, flag] = linprog(zeros(4,1), ...
    A_domain, b_modified, ones(1,4), 0, ...
    -inf(4,1), inf(4,1), options);

result.SafeDomainFeasible = flag > 0;

if flag <= 0
    return
end

%%\
%[text] ## Az új domain LP\-s ellenőrzése
%[text] A validált domainben ismét megkeressük minden CBCO pozitív és negatív irányú maximumát.
%[text]{"align":"center"} $A\_{domain}NP\\leq b\_{modified}=b\_{domain}\-IVA$
%[text] A new\_NP\_plus és new\_NP\_minus a végső maximumokat előidéző piaci pontokat tárolja.
%[text] A new\_overload az új irányonkénti maximumok FRM\-határsértése. A nulla érték azt jelenti, hogy az adott irányban nincs túllépés.
%[text] A new\_table összehasonlítja az eredeti és végső maximumokat, és melléjük rendeli az IVA\-t, valamint a végső piaci pontokat.
%[text] A SafeAfter irányonként mutatja a megfelelést. A SafeDomainVerified az all segítségével azt jelzi, hogy minden vizsgált irány megfelel\-e.
%[text] A CorrectedTable frissítése után a közvetlenül módosított irányok végső határsértése és megfelelése is látható.
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

% A kezdeti állapot, a közvetlen korrekció és a végeredmény együtt.
new_table.InitiallyUnsafe = isOverloaded;
new_table.IVA = IVA;
new_table.IVA_applied = IVA_applied;
new_table.SafeAfter = new_overload <= tol;

result.IVATable.NewFRMLevelExcess = new_overload;
result.IVATable.SafeAfter = new_overload <= tol;

result.CorrectedTable = result.IVATable(IVA_applied,:);

result.SafeResultsTable = new_table;
result.SafeOverloadedTable = new_table(new_overload > tol,:);
result.SafeDomainVerified = all(new_overload <= tol);
result.SafeLPFlags = LP_flags;

%%\
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

%[text] ## A látható korrekciók kiválasztása
%[text] A cuts azoknak a szigorított korlátoknak az indexeit tárolja, amelyek az adott metszet régi vagy új határán élként megjelennek.
%[text] A kód először ellenőrzi a jobb oldal csökkenését, majd megszámolja, hány csúcs fekszik a megfelelő határoló egyenesen.
%[text] Az nnz az igaz logikai értékeket számolja meg. Legalább két egyenesre illeszkedő csúcs jelzi a határoló élet.
%[text] A cuts rajzolási lista; az összes közvetlen IVA\-alkalmazás a result.CorrectedTable táblában található.
%[text] A plot\_rows = 0 az áttekintő ábrát jelenti. Az A–B metszetnél a cuts indexeit is hozzáadjuk, ezért ott külön részletes ábrák is készülnek.

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

%[text] ## A határoló egyenesek és az IVA
%[text] A PTDF\_row itt az adott korlát kétváltozós együtthatósora, a limit pedig a metszet új jobb oldala.
%[text] Az IVA a jobb oldalt csökkenti, miközben az együtthatósor változatlan. A régi és új egyenes ezért párhuzamos.
%[text]{"align":"center"} $b\_{2,new,i}=b\_{2,old,i}\-IVA\_i$
%[text] Ha a második együttható nem nulla, az egyenletet y\-ra rendezzük. Ellenkező esetben függőleges egyenest rajzolunk, külön kiszámított x koordinátával.

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

%[text] ## Új metszéspont és eltolási nyíl
%[text] A v1 és v2 a régi sokszög szomszédos csúcsai. A d1 és d2 előjele megmutatja, hogy a végpontok az új korlát melyik oldalán vannak.
%[text] Ha az előjelek különböznek, az él keresztezi az új határt. A kód úgy rendezi a végpontokat, hogy v1 az adott korlát megengedett oldalán legyen.
%[text] Az él pontjai:
%[text]{"align":"center"} $w=v1\+t(v2\-v1)$
%[text] A határoló egyenesbe helyettesítve számítjuk a t értékét, majd ebből az új w metszéspontot.
%[text]{"align":"center"} $t=\\frac{limit\-PTDF\_{row}v1^T}{PTDF\_{row}(v2\-v1)^T}$
%[text] A start a régi egyenesen, a w az új egyenesen van. A quiver a kettő közötti merőleges eltolást rajzolja.
%[text] A geometriai eltolás hossza az IVA és az együtthatósor normájának hányadosa.
%[text]{"align":"center"} $distance=\\frac{IVA\_i}{\\sqrt{PTDF\_{row}(1)^2\+PTDF\_{row}(2)^2}}$

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

%%\

%[appendix]
%---
