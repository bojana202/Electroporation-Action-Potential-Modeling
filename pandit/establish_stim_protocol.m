function [start_times, stim_starts, stim_ends, end_times, intervals, Istim] =...
    establish_stim_protocol(PCL, stimdelay, stimdur, stim_amp, n_stimuli, starttime_overall)

%This function takes the given stimulation parameters and creates arrays of
%time points that are used to run the simulation itself.

%start_times: the time points when every interval (equal to the duration of the
%pacing cycle length, PCL) begins. The first entry in start_times is the
%start of the simulation

%stim_starts: the time points where the injection of applied current to
%stimulate an AP begins

%stim_ends: the time points where the injection of applied current to
%stimulate an AP ends

%end_times: the time points when every interval (equal to the duration of the
%pacing cycle length, PCL) ends.  The last entry in end_times is the end of
%the simulation.

%intervals: the ode solver is called in distinct intervals.  There are
%three intervals for each cycle.  The first is from the beginning of the
%interval to the start of current injection (start_times to stim_starts).
%The second is the time of the current injection (stim_starts to
%stim_ends).  The third is the time for that cycle after the current
%injection (stim_ends ts end_times). This is a cell array that has length 
%3x of the number of cycles.  Each cell array contains the times points 
%when that interval begins and ends.  Breaking the cycle into intervals is
%necessary to ensure that the injected stimulus current begins and ends
%exactly as specified. This is 

%Istim: This is the same length as intervals, and specifies how much
%current is injected during the corresponding interval.  This will be zero
%for intervals where current is not injected, and stim_amp for those where
%it will be.

    start_times = starttime_overall + PCL*(0:n_stimuli-1);
    stim_starts = start_times + stimdelay;
    stim_ends = stim_starts + stimdur ;
    end_times = start_times + PCL ;
    
    simints = 3*n_stimuli;
    intervals = cell(1,simints);
    for i = 1:n_stimuli
        intervals{3*(i) - 2} = [start_times(i), stim_starts(i)];
        intervals{3*(i) - 1} = [stim_starts(i), stim_ends(i)];
        if i == n_stimuli
            intervals{3*i} = [stim_ends(i), end_times(end)];
        else
            intervals{3*i} = [stim_ends(i), start_times(i+1)];
        end
    end
    Istim = zeros(simints,1) ;
    stimindices = 3*(1:n_stimuli) - 1 ;
    Istim(stimindices) = -stim_amp ;

return
