function [redispatch_result_table] = redispatch(cbco_table, ref_exp, ...
    gsk_for_lsqlin, n2sPTDF_cont, n2sPTDF_for_lsqlin, presolved_cbco_table, zones)


    % minden csomóponthoz hozzárendelünk egy deltát (döntési változó)
    
    % első eset: bármelyik csomópont bármennyivel
    
    % np = node2zone * x
    node2zone = (gsk_for_lsqlin>0)';
    
    % F0 = Fref - z2sPTDF * NP
    % F0 = n2sPTDF*ref_exp - z2sPTDF*node2zone*ref_exp
    % F0 = (n2sPTDF - z2sPTDF*node2zone)*ref_exp
    % F0(delta) = (n2sPTDF - z2sPTDF*node2zone)*(ref_exp - delta)
    
    % ne legyen amr
    % RAM0 = Fmax-FRM-F0 >= minRAM
    % F0(delta) <= Fmax - FRM - minRAM
    % (n2sPTDF - z2sPTDF*node2zone)*(ref_exp - delta) <= Fmax - FRM - minRAM
    % (z2sPTDF*node2zone - n2sPTDF) * delta <= 
    % Fmax - FRM - minRAM - (n2sPTDF - z2sPTDF*node2zone)*ref_exp
    
    z2sPTDF_redispatch = table2array(presolved_cbco_table(:,3:6));
    
    % n2sPTDF a teljes cbco-ra vonatkozik
    n2sPTDF_full = [n2sPTDF_for_lsqlin; -n2sPTDF_for_lsqlin; n2sPTDF_cont];
    
    [~,id] = ismember(presolved_cbco_table, cbco_table, "rows");
    n2sPTDF_redispatch = n2sPTDF_full(id, :);

    A = z2sPTDF_redispatch * node2zone - n2sPTDF_redispatch;
    A_lp = [A, -A];

    Fmax   = presolved_cbco_table.Fmax;
    FRM    = presolved_cbco_table.FRM;
    minRAM = presolved_cbco_table.minRAM;
    
    b = Fmax - FRM - minRAM - (n2sPTDF_redispatch - z2sPTDF_redispatch * node2zone) * ref_exp;
    
    % delták összege 1 (negatívok is pozitív értéket vesznek fel, emiatt -)
    Aeq = [ones(1, height(ref_exp)), -ones(1, height(ref_exp))];
    beq = 0;
    
    % delta+, delta- csak pozitív értéket kaphat
    lb_case1 = zeros(2*height(ref_exp),1);
    ub_case1 = inf(2*height(ref_exp),1);
    f = ones(2*height(ref_exp),1);

    options = optimoptions("linprog", "Display", "none");
    x_case1 = linprog(f,A_lp,b,Aeq,beq,lb_case1,ub_case1,options);
    
    % redispatch_case1 = x_case1;

    redispatch_case1 = x_case1(1:height(ref_exp)) - x_case1(height(ref_exp)+1:end);
    

    % második eset: bármelyik csomópont korlátozott deltával/pmin és pmax
    % deltára alsó és felső korlát
    % lehetne pmin-re: p-delta-
    % lehetne pmax-ra: p+delta+
    
    lb_case2 = zeros(2*height(ref_exp),1);
    ub_case2 = 35 * ones(2*height(ref_exp),1);
 
    x_case2 = linprog(f,A_lp,b,Aeq,beq,lb_case2,ub_case2,options);

    % redispatch_case2 = x_case2;
    redispatch_case2 = x_case2(1:height(ref_exp)) - x_case2(height(ref_exp)+1:end);
    
    % harmadik eset: csak adott csomópontok vehetnek részt:
    % a nem résztvevő csomópontok deltája 0-ra rögzítve:
    % N1, N2, N3, N7, N9, N11, N14 nem szolgáltat redispatch-et
    
    lb_case3 = lb_case2;
    ub_case3 = ub_case2;
    ub_case3([1 2 3 7 9 11 14]) = 0;
    ub_case3(height(ref_exp) + [1 2 3 7 9 11 14]) = 0;

    x_case3 = linprog(f,A_lp,b,Aeq,beq,lb_case3,ub_case3,options);

    % redispatch_case3 = x_case3;
    redispatch_case3 = x_case3(1:height(ref_exp)) - x_case3(height(ref_exp)+1:end);

    % negyedik eset: csak zónabeli csomópontok között mehet végbe (akár több zónában egymás
    % mellett párhuzamosan)
    % deltákra zónás egyensúlyi egyenletek

    lb_case4 = lb_case2;
    ub_case4 = ub_case2;

    Aeq_case4 = [node2zone, -node2zone];
    beq_case4 = zeros(numel(zones),1);

    x_case4 = linprog(f,A_lp,b,Aeq_case4,beq_case4,lb_case4,ub_case4,options);

    % redispatch_case4 = x_case4;
    redispatch_case4 = x_case4(1:height(ref_exp)) - x_case4(height(ref_exp)+1:end);

    redispatch_result_table = table(redispatch_case1, redispatch_case2, redispatch_case3, redispatch_case4, 'VariableNames',{'Case1', 'Case2', 'Case3', 'Case4'});

end