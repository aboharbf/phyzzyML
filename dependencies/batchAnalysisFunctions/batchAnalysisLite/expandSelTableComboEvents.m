function selTable = expandSelTableComboEvents(selTable, params);
% Inputs:
% selTable - usual output of selectivity functions in runAnalyses
% params - struct containing:
% .comboEvents - the name of the combo event to be created n*1
% .comboSubEvents - a cell array, where each entry is a cell array of the
% events represented in selTable to be combined to make the combo event.

% Create Combo event columns in the selTable
comboEvents = params.comboEvents;
comboSubEvents = params.comboSubEvents;

varNames = selTable.Properties.VariableNames;
for combo_i = 1:length(comboEvents)
  % Initialize combo index, and collect subEvents.
  subEvents = comboSubEvents{combo_i};
  comboInd = false(size(selTable,1),1);
  
  % See if individual events are present
  for sub_i = 1:length(subEvents)
    if any(strcmp(varNames, subEvents{sub_i}))
      comboInd = comboInd | (selTable.(subEvents{sub_i}) ~= 0);
    end
  end
  
  % Add to larger table
  selTable.(comboEvents{combo_i}) = comboInd;
  
end