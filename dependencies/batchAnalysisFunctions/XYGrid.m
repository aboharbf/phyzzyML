function figHand = XYGrid(dataTable, RowLabels, ColumnLabels, params)
% Function which generates a figure with a grid representing non-empty
% entries in a table or cell array (or 0's in a numeric array).
% Inputs:
% dataTable - a cell array of values where the rows are the X labels, and
% the columns are the Y labels. Can also just be a table.
% X Labels - a vector of names matching the first dimension of the dataTable.
% Y Labels - a vector of names matching the second dimension of the
% dataTable.
% params - a struct containing the fields:

if isa(dataTable, 'table')
  % If input data is a table, extract needed info.
  
  if isempty(RowLabels)
    RowLabels = dataTable.Properties.RowNames;
  end
  
  if isempty(ColumnLabels)
    ColumnLabels = dataTable.Properties.VariablesNames;
  end
  
  dataMat = cellfun(@(x) size(x, 1), table2cell(dataTable));
  
else
  dataMat = dataTable;
end


figHand = figure('Units', 'Normalized');
axHand = axes();
imagesc(dataMat);

% The main problem this code needs to solve is updating the values XTick,
% YTick, XLabel, and YLabel according to zoom.

% Load the full Labels into the user data of the figure.
axHand.UserData.ColumnLabels = ColumnLabels;
axHand.UserData.RowLabels = RowLabels;
axHand.TickLength = [0.003 0.03];
axHand.XTickLabelRotation = 45;

axHand.XTick = 1:ceil(length(ColumnLabels)/30):length(ColumnLabels);
axHand.YTick = 1:ceil(length(RowLabels)/30):length(RowLabels);

axHand.XTickLabel = ColumnLabels(axHand.XTick);
axHand.YTickLabel = RowLabels(axHand.YTick);

% Add Colorbar
colorbar()

% Add the needed callback handles to appropriately shift the axes.
[figHand.WindowButtonUpFcn, figHand.SizeChangedFcn] = deal({@adjustAxesLabels});

% Add the labels, if present.
if isfield(params, 'XLabel')
  xlabel(params.XLabel);
end

if isfield(params, 'YLabel')
  ylabel(params.YLabel);
end

if isfield(params, 'plotTitle')
  title(params.plotTitle);
end

end

function adjustAxesLabels(src,event)
  % Callback for when function is zoomed. Uses structures in the user data
  % to adjust X and Y Tick and TickLabel.
  
  LabelDensity = 40; % Number of Labels per normalized unit of size in the plot.
  axHand = findobj(src.Children, 'Type', 'Axes');
  
  % Use a combination of the normalized units the figure takes up + the
  % number of labels in the current zooming of the figure to decide how
  % many to show.
  axX = [round(axHand.XLim(1)) floor(axHand.XLim(2))];
  axY = [round(axHand.YLim(1)) floor(axHand.YLim(2))];
  
  axSize = axHand.Position(3:4);
  axRange = [diff(axX) diff(axY)];
  axDens = axRange ./ (LabelDensity * axSize);
  axStep = round(max(axDens, 1));
    
  axHand.XTick = axX(1):axStep(1):axX(2);
  axHand.YTick = axY(1):axStep(2):axY(2);

  axHand.XTickLabel = axHand.UserData.ColumnLabels(axHand.XTick');
  axHand.YTickLabel = axHand.UserData.RowLabels(axHand.YTick);

end