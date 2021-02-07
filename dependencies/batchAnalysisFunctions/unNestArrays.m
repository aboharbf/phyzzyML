function OutCellArray = unNestArrays(inCell)
% Function which unnested/nests a cell array into a single huge matrix, and
% stores the borders between the original cells.



% First Input can either be
% - an array, where each number represents the max index at that tier
% OR
% - a nested cell array, which will be copied in structure
% 2nd Argument = 'ones','zeros','NaN', or 'cell', to distinguish structure
% at bottom. default is cell
% 3rd argument = array of dimensions to feed into 2nd argument function.
% default is [0, 0].
% 4th argument = depth, if the full depth isn't desired. (default = 100)

%Defaults: (Array, 'cell', [0 0], 100)

%% Function
if ~iscell(inCell{1})
  % What to do once you've reached the bottom.
  OutCellArray = vertcat(inCell{:});
else
  OutCellArray = cell(size(inCell));
  for inCell_ind = 1:length(inCell)
    OutCellArray{inCell_ind} = unNestArrays(inCell{inCell_ind});
  end
  OutCellArray = vertcat(OutCellArray{:});
end