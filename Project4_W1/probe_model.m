%% probe_model.m — inspect Battery block connections

load_system('p4_model.slx');

batPath = 'p4_model/Longitudinal Dynamics Body Frame/Subsystem/Electric Motor Drive/Battery';

% Get Unit Delay details
udPath = [batPath '/Unit Delay'];
fprintf('=== Unit Delay ===\n');
fprintf('  InitialCondition: %s\n', get_param(udPath, 'InitialCondition'));
ph = get_param(udPath, 'PortHandles');
fprintf('  Inport connections:\n');
for i = 1:length(ph.Inport)
    line = get_param(ph.Inport(i), 'Line');
    if line > 0
        srcPort = get_param(line, 'SrcPortHandle');
        srcBlock = get_param(get_param(srcPort,'Parent'), 'Name');
        fprintf('    <- %s\n', srcBlock);
    end
end
fprintf('  Outport connections:\n');
for i = 1:length(ph.Outport)
    line = get_param(ph.Outport(i), 'Line');
    if line > 0
        dstPorts = get_param(line, 'DstPortHandle');
        for j = 1:length(dstPorts)
            dstBlock = get_param(get_param(dstPorts(j),'Parent'), 'Name');
            fprintf('    -> %s\n', dstBlock);
        end
    end
end

% Trace the Sum (algebraic variable) block
fprintf('\n=== Sum block ===\n');
sumPath = [batPath '/Sum'];
ph2 = get_param(sumPath, 'PortHandles');
fprintf('  Inputs:\n');
for i = 1:length(ph2.Inport)
    line = get_param(ph2.Inport(i), 'Line');
    if line > 0
        srcPort = get_param(line, 'SrcPortHandle');
        srcBlock = get_param(get_param(srcPort,'Parent'), 'Name');
        fprintf('    [port %d] <- %s\n', i, srcBlock);
    end
end
fprintf('  Output:\n');
for i = 1:length(ph2.Outport)
    line = get_param(ph2.Outport(i), 'Line');
    if line > 0
        dstPorts = get_param(line, 'DstPortHandle');
        for j = 1:length(dstPorts)
            dstBlock = get_param(get_param(dstPorts(j),'Parent'), 'Name');
            fprintf('    -> %s\n', dstBlock);
        end
    end
end

% Trace I_battery inport
fprintf('\n=== I_battery inport ===\n');
iPath = [batPath '/I_battery'];
ph3 = get_param(iPath, 'PortHandles');
for i = 1:length(ph3.Outport)
    line = get_param(ph3.Outport(i), 'Line');
    if line > 0
        dstPorts = get_param(line, 'DstPortHandle');
        for j = 1:length(dstPorts)
            dstBlock = get_param(get_param(dstPorts(j),'Parent'), 'Name');
            fprintf('  -> %s\n', dstBlock);
        end
    end
end

% Check Battery subsystem atomic setting
fprintf('\n=== Battery subsystem atomic params ===\n');
try; fprintf('  TreatAsAtomicUnit: %s\n', get_param(batPath,'TreatAsAtomicUnit')); catch; end
try; fprintf('  MinAlgLoopOccurrences: %s\n', get_param(batPath,'MinAlgLoopOccurrences')); catch; end

% Get algebraic loops
fprintf('\n=== Algebraic loops ===\n');
loops = Simulink.BlockDiagram.getAlgebraicLoops('p4_model');
for i = 1:length(loops)
    fprintf('Loop %d blocks:\n', i);
    for j = 1:length(loops(i).Blocks)
        fprintf('  %s\n', loops(i).Blocks{j});
    end
end
