function data = removeDups(data)
% for 3 runs on 06/23/2021, I left in duplicates of eventmarkers put there
% for the sake of testing a problem with eventmarker transmission to
% blackrock. Those weren't removed prior to the runs, so now these 3 runs
% have duplicates. This functions removes them and passes the data back. 

dupEvents = [10 20 30 40];

for trial_i = 1:length(data)
  codes2Update = data(trial_i).BehavioralCodes;
  keepInd = false(size(codes2Update.CodeNumbers));
  
  % Keep the first 2 and last markers, and any errors.
  keepInd(1:2) = true;
  keepInd(end) = true;
  keepInd(codes2Update.CodeNumbers == 3) = true;
  keepInd(codes2Update.CodeNumbers == 4) = true;
  
  % Find the first of the duplicatedEvents, and keep that one
  for ev_i = 1:length(dupEvents)
    keepInd(find(codes2Update.CodeNumbers == dupEvents(ev_i), 1)) = true;
  end
  
  % Update the structures, pass back
  codes2Update.CodeNumbers = codes2Update.CodeNumbers(keepInd);
  codes2Update.CodeTimes = codes2Update.CodeTimes(keepInd);
  data(trial_i).BehavioralCodes = codes2Update;
  
end


