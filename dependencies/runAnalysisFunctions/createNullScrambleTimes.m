function nullEventTimes = createNullScrambleTimes(eventTimes, maxShift, minEventDist)
% Takes a list of event times, and generates a list of scramble times based
% on input parameters, in a way which isn't too close to any other events.
% Inputs:
% - eventTimes: a cell array for every event which needs to be scrambled.

nullEventTimes = cell(size(eventTimes));
for event_i = 1:length(eventTimes)
  
  % Extract information, find the earliest and latest times
  eventsTiled = eventTimes{event_i};
  firstTime = min(eventsTiled);
  lastTime = max(eventsTiled);

  % Generate a random vector of shifts.
  shiftCount = length(eventsTiled)*10;
  negShifts = randi([-maxShift, -minEventDist], [shiftCount/2, 1]);
  posShifts = randi([minEventDist, maxShift], [shiftCount/2, 1]);
  shifts = [negShifts; posShifts];
  shifts = shifts(randperm(length(shifts)));
  
  % Add it to an equally sized mix of the events themselves.
  nullTimesForEvent = eventsTiled(randi(length(eventsTiled), [shiftCount,1])) + shifts;
  
  % check if any of these null times are too close to event times.
  nullTimesForEventCheck = repmat(nullTimesForEvent, [1, length(eventsTiled)]);
  nullTimesForEventCheck = (nullTimesForEventCheck' - eventsTiled)';
  nullTimeTooCloseInd = abs(nullTimesForEventCheck) <= minEventDist;
  keepInd = ~any(nullTimeTooCloseInd,2);
  nullTimesForEvent = nullTimesForEvent(keepInd);
  
  % Remove ones that are after the last time stamp or before the first
  keepInd = nullTimesForEvent > firstTime & nullTimesForEvent < lastTime;
  nullTimesForEvent = nullTimesForEvent(keepInd);
  
  % Store into appropriate cell array
  nullEventTimes{event_i} = nullTimesForEvent;
end

end