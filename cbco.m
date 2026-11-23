function [cbco_table] = cbco(z2sPTDF, LODF, G, Edge_table, f_max, min_ram, frm, zones, ref_net_pos)

    % alapeset
    z2sPTDF_basecase = [z2sPTDF; -z2sPTDF];

    [~,loc] = ismember(Edge_table.Name, G.Edges.Name);
    cb_basecase = repmat(G.Edges.Name(loc),2,1);
    co_basecase = repmat("BaseCase", height(z2sPTDF_basecase),1);
    cbco_table_basecase = table(cb_basecase, co_basecase, 'VariableNames', {'CB', 'CO'});
    
    for z = 1:numel(zones)
        col_name = "z2sPTDF_" + zones(z);
        cbco_table_basecase.(col_name) = z2sPTDF_basecase(:, z);
    end
    f_ref_basecase = G.Edges.("Reference flow")(loc);
    f_0_basecase = G.Edges.("Zero balance reference flow")(loc);
    ram0_p_basecase = f_max - frm - f_0_basecase;
    ram0_m_basecase = f_max - frm + f_0_basecase;
    amr_p_basecase = max(0, min_ram - ram0_p_basecase);
    amr_m_basecase = max(0, min_ram - ram0_m_basecase);
    ram_p_basecase = ram0_p_basecase + amr_p_basecase;
    ram_m_basecase = ram0_m_basecase + amr_m_basecase;
  
    cbco_table_basecase.RAM0 = [ram0_p_basecase; ram0_m_basecase];
    cbco_table_basecase.AMR = [amr_p_basecase; amr_m_basecase];
    cbco_table_basecase.RAM = [ram_p_basecase; ram_m_basecase];

    % kontingenciák

    % összes lehetséges CBCO kombináció összeállítása
    cbco_combs = combinations(G.Edges.Name(loc), G.Edges.Name(loc));
    cbco_combs.Properties.VariableNames = {'CO','CB'};
    cbco_combs(cbco_combs.CB == cbco_combs.CO,:) = [];
    [~,cb_id] = ismember(cbco_combs.CB, G.Edges.Name(loc)); 
    [~,co_id] = ismember(cbco_combs.CO, G.Edges.Name(loc));

    % LODF értékek kiválasztása minden CBCO párhoz
    lodf_cbco = zeros(numel(cb_id), 1);
    for i = 1:numel(cb_id)
        lodf_cbco(i) = LODF(cb_id(i), co_id(i));
    end

    % elemenkénti szorzással
    % pl.: L1 kiesése
    % z2sPTDF_L1 = z2sPTDF + LODF(:,1) .* z2sPTDF(1,:);
    % z2sPTDF_L1(1,:) = []
   
    z2sPTDF_cont = z2sPTDF(cb_id,:) + lodf_cbco .* z2sPTDF(co_id,:);

    f_ref_cont = f_ref_basecase(cb_id) + lodf_cbco .* f_ref_basecase(co_id);
    f_0_cont = f_ref_cont - z2sPTDF_cont * ref_net_pos(:);
    
    % pozitív és negatív irányokhoz külön sor
    cb_cbco = [cbco_combs.CB; cbco_combs.CB];
    co_cbco = [cbco_combs.CO; cbco_combs.CO];
    
    cbco_table_cont = table(cb_cbco, co_cbco, 'VariableNames', {'CB', 'CO'});
    
    % z2sPTDF-ek hozzáadása
    z2sPTDF_cont = [z2sPTDF_cont; -z2sPTDF_cont];
    for z = 1:numel(zones)
        col_name = "z2sPTDF_" + zones(z);
        cbco_table_cont.(col_name) = z2sPTDF_cont(:,z);
    end
    
    ram0_p_cont = f_max - frm - f_0_cont;
    ram0_m_cont = f_max - frm + f_0_cont;
    
    amr_p_cont = max(0, min_ram - ram0_p_cont);
    amr_m_cont = max(0, min_ram - ram0_m_cont);
    
    ram_p_cont = ram0_p_cont + amr_p_cont;
    ram_m_cont = ram0_m_cont + amr_m_cont;

    cbco_table_cont.RAM0 = [ram0_p_cont; ram0_m_cont];
    cbco_table_cont.AMR = [amr_p_cont; amr_m_cont];
    cbco_table_cont.RAM = [ram_p_cont; ram_m_cont];
    
    cbco_table = [cbco_table_basecase; cbco_table_cont];

    % fmax, minram, fref, f0
    cbco_table.Fmax = f_max * ones(height(cbco_table),1);
    cbco_table.minRAM = min_ram * ones(height(cbco_table),1);

    % negatív irányban (-1)-gyel szorozva a flow, igazítva ezzel a
    % korlátok kényelmes felírásához
    cbco_table.Fref = [f_ref_basecase; -f_ref_basecase; f_ref_cont; -f_ref_cont];
    cbco_table.F0 = [f_0_basecase; -f_0_basecase; f_0_cont; -f_0_cont]
end