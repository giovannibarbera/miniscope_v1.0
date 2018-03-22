function [CA_TRACES BEHAVIOR] = load_mouse(wspace, cell_type, mouse_numb, N_days, N_exp, N)
%LOAD_MOUSE Load all data relative to a specific mouse, including behavior
%
% cell_type:        (1), 1 for D1 neurons, 2 for D2 neurons
% mouse_numb:       (1), mouse number
% day:              (1), day number
% exper:            (1), experiment number
% N_days:           (1), number of days
% N_exp:            (1), number of experiments
% N:                (1), number of datapoints per experiment (time length)
%
% CA_TRACES:        (N_days, N_exp, N, N_cells), raw calcium traces
% BEHAVIOR:         (cell) behavior vector as defined in BehaviorDataRead.m
% 
%   giobarbera@neuralmapper.com



eval(['tot_cells = load(wspace, ''identified_cells_daily_all_',num2str(cell_type),'_',num2str(mouse_numb),''');'])
eval(['tot_cells = tot_cells.identified_cells_daily_all_',num2str(cell_type),'_',num2str(mouse_numb),';'])  
eval(['BEHAVIOR = load(wspace, ''behavior_',num2str(cell_type),'_',num2str(mouse_numb),''');'])
eval(['BEHAVIOR = BEHAVIOR.behavior_',num2str(cell_type),'_',num2str(mouse_numb),';'])
eval(['calcium_traces = load(wspace, ''firing_cells_daily_all_',num2str(cell_type),num2str(mouse_numb,'%02i'),'*'');'])


CA_TRACES = zeros(N_days, N_exp, N, tot_cells);


for curr_day = 1:N_days
    for curr_exper = 1:N_exp
        if curr_exper<=3 
            aw_co = 'A';
        else 
            aw_co = 'C';
        end

        curr_header = strcat(num2str(cell_type),num2str(mouse_numb,'%02i'), num2str(curr_day,'%02i'),num2str(curr_exper,'%02i'),aw_co);


        % Calcium traces
        eval(['CA_TRACES(curr_day,curr_exper,:,:) = calcium_traces.firing_cells_daily_all_',curr_header,';'])
    end
end


