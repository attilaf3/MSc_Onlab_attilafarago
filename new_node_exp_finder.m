function [node_ref_exp_new] = new_node_exp_finder(cbco_table, ...
    node_ref_exp_old, n2sPTDF_for_lsqlin, gsk_for_lsqlin, n2sPTDF_cont)

    C = eye(height(node_ref_exp_old));
    d = [node_ref_exp_old];
      
    % GSK nemnulla elemei 1-re állítása
    node2zone = (gsk_for_lsqlin>0)';
    % áttérés a csomóponti vektor és az NP között mátrixszorzással
    % np = node2zone * x

    % F0 = Fref - z2sPTDF * NP
    % F0 = n2sPTDF*x - z2sPTDF*node2zone*x
    % F0 = (n2sPTDF - z2sPTDF*node2zone)*x

    % alapesetben ne legyen amr
    % tehát RAM0 = Fmax-FRM-F0 >= minRAM
    % F0 <= Fmax - FRM - minRAM
    % (n2sPTDF - z2sPTDF*node2zone)*x <= Fmax - FRM - minRAM

    % bal oldal
    % n2sPTDF_for_lsqlin csak a pozitív irányra vonatkozik alapból
    % a z2sPTDF_base tartalmazza mindkét irányt a cbco fv miatt 
    n2sPTDF_base = [n2sPTDF_for_lsqlin; -n2sPTDF_for_lsqlin];
    z2sPTDF_base = table2array(cbco_table(cbco_table.CO == "BaseCase", 3:6));

    % jobb oldal
    Fmax_base = cbco_table.Fmax(cbco_table.CO == "BaseCase");
    FRM_base = cbco_table.FRM(cbco_table.CO == "BaseCase");
    minRAM_base = cbco_table.minRAM(cbco_table.CO == "BaseCase");


    % n-1 esetben ne alakuljon ki (fmax-frm)-nél nagyobb áramlás egyik vezetéken
    % sem (amr viszont megengedett)
    % Fref <= Fmax-FRM
    % n2sPTDF * x <= Fmax -FRM
    % bal oldal
    n2sPTDF_cont;

    % jobb oldal
    Fmax_cont = cbco_table.Fmax(cbco_table.CO ~= "BaseCase");
    FRM_cont = cbco_table.FRM(cbco_table.CO ~= "BaseCase");

    % egyenlőtlenségi korlátok
    A = [n2sPTDF_base - z2sPTDF_base * node2zone; n2sPTDF_cont];
    b = [Fmax_base - FRM_base - minRAM_base; Fmax_cont - FRM_cont];
    
    % egyenlőség
    Aeq = ones(1, height(node_ref_exp_old));
    beq = 0;

    lb = -10000*ones(height(node_ref_exp_old),1);
    ub = 10000*ones(height(node_ref_exp_old),1);
    
    % x = lsqlin(C,d,A,b,Aeq,beq,lb,ub);
    % toolbox/verzió hiba miatt másik algoritmust kellett beállítani
    options = optimoptions('lsqlin', 'Algorithm', 'active-set', 'Display', 'iter');
    x0 = zeros(size(d));
    x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options);

    node_ref_exp_new = x;

end