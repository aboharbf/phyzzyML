function saveEnv(varargin)
% Function saves the environment, specifically environments with the
% spikeDataBank, which is too large to save in a normal file w/o excessive
% compression time.

if nargin == 0
  saveSpikeDB = 1;
elseif nargin == 1
  saveSpikeDB = varargin{1};
elseif nargin > 1
  error('Too many arguments');  
end


  % Parse variables of the environment, pick out non-spikeDataBank ones. 
  vars = evalin('caller','whos');
  varNames = {vars(:).name};
  varNames(strcmp('spikeDataBank',varNames)) = [];
  
  % Assign the shortened list in the environment, then save one and then
  % the other.
  assignin('caller','varNames', varNames)
  evalin('caller', "save(fullfile(outputDir, [preprocessParams.spikeDataFileName 'Vars']), varNames{:})")
  if saveSpikeDB
    evalin('caller', "saveSpikeDataBank(spikeDataBank, 2, 'save',outputDir)" )
  end
end