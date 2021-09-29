function tableArrayClean = removeTableVariable(tableArray, var2remove)
% A function which removes a table variable if present for a set of tables
% saved in a cell array.

tableArrayClean = cell(size(tableArray));
for ii = 1:length(tableArray)
  % Pull table
  tableTmp = tableArray{ii};
  var2RemoveTmp = intersect(tableTmp.Properties.VariableNames, var2remove);
  
  % Empty of undesired variables
  for var_i = 1:length(var2RemoveTmp)
    tableTmp.(var2RemoveTmp{var_i}) = [];
  end
  
  % Return table
  tableArrayClean{ii} = tableTmp;
  
end
end