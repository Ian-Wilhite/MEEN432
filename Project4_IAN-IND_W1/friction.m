function [BPP_friction, BPP_regen] = friction_brake_limit(BPP_total, v_mps)
%FRICTION_BRAKE_LIMIT  Enforce minimum friction-braking percentage (Project 4 req.)
%
%  The rules:
%    v < 5 mph  : 100% friction,  0% regen allowed
%    5–25 mph   : minimum friction drops linearly from 100% → 5%
%    v > 25 mph : minimum 5% friction, up to 95% regen allowed
%
%  Inputs
%    BPP_total  – total brake pedal position [0–1] (combined friction + regen)
%    v_mps      – current vehicle speed (m/s)
%
%  Outputs
%    BPP_friction – fraction of braking done by friction pads [0–1]
%    BPP_regen    – fraction available for regenerative braking [0–1]

    mph2mps = 0.44704;
    v_mph   = v_mps / mph2mps;

    % Minimum required friction-braking FRACTION of total brake demand
    if v_mph < 5
        min_fric_frac = 1.0;
    elseif v_mph <= 25
        % Linear interpolation: 1.0 at 5 mph → 0.05 at 25 mph
        min_fric_frac = 1.0 - 0.95 * (v_mph - 5) / 20;
    else
        min_fric_frac = 0.05;
    end

    BPP_friction = BPP_total * min_fric_frac;
    BPP_regen    = BPP_total * (1 - min_fric_frac);
end
