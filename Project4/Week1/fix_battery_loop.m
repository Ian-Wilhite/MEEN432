%% fix_battery_loop.m
% Breaks the algebraic loop in the Battery subsystem by wiring the
% existing (orphaned) Unit Delay into the Gain -> Sum feedback path.
%
% Before: Gain --> Sum port 2          (algebraic loop)
% After:  Gain --> Unit Delay --> Sum port 2   (loop broken)

load_system('p4_model.slx');

batPath = 'p4_model/Longitudinal Dynamics Body Frame/Subsystem/Electric Motor Drive/Battery';
gainPath = [batPath '/Gain'];
sumPath  = [batPath '/Sum'];
udPath   = [batPath '/Unit Delay'];

fprintf('Checking existing connections...\n');

% Verify Unit Delay is truly unconnected
udPH = get_param(udPath, 'PortHandles');
ud_in_line  = get_param(udPH.Inport(1),  'Line');
ud_out_line = get_param(udPH.Outport(1), 'Line');
assert(ud_in_line  < 0, 'Unit Delay inport already connected — check model manually.');
assert(ud_out_line < 0, 'Unit Delay outport already connected — check model manually.');
fprintf('  Unit Delay confirmed unconnected.\n');

% Find the line from Gain output -> Sum port 2
gainPH = get_param(gainPath, 'PortHandles');
gain_out_line = get_param(gainPH.Outport(1), 'Line');
assert(gain_out_line > 0, 'Gain output has no line — unexpected.');

% Confirm it goes to Sum port 2
dstPorts = get_param(gain_out_line, 'DstPortHandle');
assert(length(dstPorts) == 1, 'Gain output fans out to multiple destinations — fix manually.');
dstBlock = get_param(get_param(dstPorts(1), 'Parent'), 'Name');
assert(strcmp(dstBlock, 'Sum'), 'Gain does not connect to Sum — topology changed, fix manually.');
fprintf('  Confirmed: Gain -> Sum (port 2) line exists.\n');

% Get Sum port 2 handle
sumPH = get_param(sumPath, 'PortHandles');
sum_port2 = sumPH.Inport(2);

% Position Unit Delay between Gain and Sum
gain_pos = get_param(gainPath, 'Position');   % [left top right bottom]
sum_pos  = get_param(sumPath,  'Position');
ud_w = 30; ud_h = 30;
ud_cx = round((gain_pos(3) + sum_pos(1)) / 2);
ud_cy = round((gain_pos(2) + gain_pos(4)) / 2);
new_ud_pos = [ud_cx - ud_w/2, ud_cy - ud_h/2, ud_cx + ud_w/2, ud_cy + ud_h/2];
set_param(udPath, 'Position', new_ud_pos);
fprintf('  Repositioned Unit Delay to [%d %d %d %d].\n', new_ud_pos);

% Delete the old Gain -> Sum line
delete_line('p4_model/Longitudinal Dynamics Body Frame/Subsystem/Electric Motor Drive/Battery', ...
    'Gain/1', 'Sum/2');
fprintf('  Deleted old Gain -> Sum line.\n');

% Add Gain -> Unit Delay
add_line('p4_model/Longitudinal Dynamics Body Frame/Subsystem/Electric Motor Drive/Battery', ...
    'Gain/1', 'Unit Delay/1', 'autorouting', 'on');
fprintf('  Added Gain -> Unit Delay.\n');

% Add Unit Delay -> Sum port 2
add_line('p4_model/Longitudinal Dynamics Body Frame/Subsystem/Electric Motor Drive/Battery', ...
    'Unit Delay/1', 'Sum/2', 'autorouting', 'on');
fprintf('  Added Unit Delay -> Sum.\n');

% Save the model
save_system('p4_model.slx');
fprintf('\nModel saved. Algebraic loop fix applied.\n');
fprintf('Verify by running p4_runsim — the algebraic loop warning should be gone.\n');
