% betöltöm a cbco_table-t
% inicializálok egy presolved_cbco_table-t ugyanazokkal az oszlopnevekkel
% végig iterálok az összes cbco-n
% az adott CBCO-t úgy vizsgálom, hogy a baloldala lesz az adott LP
% célfüggvénye (maximalizálás), míg a korlátok az összes többi CBCO + NP
% egyensúly
% ezt követően a célfüggvény értéket összevetem a vizsgált CBCO jobb
% oldalával
% ha nagyobb akkor a célfüggvényérték, akkor nem redundáns, ekkor mehet a
% presolved_cbco_table új sorába

function presolved_cbco_table = presolver(cbco_table)
   presolved_cbco_table = table('VariableNames',cbco_table.Properties.VariableNames)

end

%[appendix]
%---
