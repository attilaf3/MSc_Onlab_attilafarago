function [cbco_with_iva] = validation(cbco_table, zones)

    cbco_with_iva = cbco_table;
    iva = zeros(height(cbco_with_iva),1);
    % f0 = fref-z2sPTDF*NP
    % adott NP-kre:
    % f = f0 + z2sptdf*np
    % végigmegyünk az összes cbco-n (presolved)
    % és maximalizáljuk az adott f(cbco) = f0+z2sPTDF*NP(CBCO)-t
    % betartva a korlátokat
    % majd pedig a túlterhelést nézzük meg, és a legnagyobbat összevetjük a
    % RAM0-val
    overload = zeros(height(cbco_with_iva),1);
    cbco_with_iva.IVA = iva;
    A = table2array(cbco_with_iva(:,3:6));

    Aeq = ones(1,numel(zones));
    beq = 0;

    while true  
        
        % jobb oldal
        % ram0 + amr - iva = ram - iva
        b = cbco_with_iva.RAM - cbco_with_iva.IVA;

        for row=1:height(cbco_with_iva)
    
            f = table2array(cbco_with_iva(row,3:6));
       
            options = optimoptions('linprog', 'Display', 'none');
            [~, fval, ~] = linprog(-f, A, b, Aeq, beq, [], [], options);
        
            overload(row) = - fval - cbco_with_iva.RAM0(row);
        end
        
        max_overload = max(overload);
        max_overload_id = find(overload == max_overload);
        
        % 0.001, különben beragadna a ciklus
        if max_overload > 0.001
            cbco_with_iva.IVA(max_overload_id) = cbco_with_iva.IVA(max_overload_id) + max_overload;
        else
            break
        end
    end

end