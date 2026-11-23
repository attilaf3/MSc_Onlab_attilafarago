function [presolved_cbco_table] = presolve(cbco_table, zones)

% korlátok
% ezek egyben
% z2sPTDF * NP <= RAM+
% -z2sPTDF * NP <= RAM-
% sumNP = 0

Aeq = ones(1,numel(zones));
beq = 0;

% presolved_cbco_ba bekerülő cbco-k indexe
not_redundant_rows = [];

% iterálás közben a redundánsnak talált sorok elhagyása
redundant_rows = [];

% végig iterálok az összes cbco-n
    for cbco_id = 1:height(cbco_table)
        z2sPTDF_cbco = [cbco_table.z2sPTDF_A cbco_table.z2sPTDF_B cbco_table.z2sPTDF_C cbco_table.z2sPTDF_D]; 
        RAM_cbco = cbco_table.RAM;

        % az adott CBCO-t úgy vizsgálom, hogy a baloldala lesz az adott LP
        vizsgalt_a = z2sPTDF_cbco(cbco_id,:);
        vizsgalt_b = RAM_cbco(cbco_id);    
  
        % célfüggvénye (maximalizálás), míg a korlátok az összes többi CBCO + NP
        z2sPTDF_cbco([redundant_rows, cbco_id],:) = [];
        RAM_cbco([redundant_rows, cbco_id]) = [];

        A = [z2sPTDF_cbco];
        b = [RAM_cbco];
        f = [-vizsgalt_a];
        [~, fval, ~] = linprog(f, A, b, Aeq, beq);

        if -fval >= vizsgalt_b
            not_redundant_rows(end+1) = cbco_id;
        else
            redundant_rows(end+1) = cbco_id;
        end
    end
presolved_cbco_table = cbco_table(not_redundant_rows, :);

end