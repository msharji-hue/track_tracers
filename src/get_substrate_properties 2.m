function sub = get_substrate_properties(material, condition)
% GET_SUBSTRATE_PROPERTIES  Single source of truth for measured substrate bed
% properties, keyed on (material, condition). Mirrors get_calibration.m.
%
%   sub = get_substrate_properties('GB','full')
%
%   Returns: material, condition, rho_particle_g_cm3, rho_bulk_g_cm3, phi,
%            ok, reason.
%
%   Definitions / usage:
%     rho_particle : true (skeletal) particle density — the phi reference only.
%                    (soda-lime glass 2.50; Hess pumice / amorphous aluminum
%                    silicate 2.35 from SDS specific gravity)
%     rho_bulk     : measured bed bulk density (bed mass / bed volume). THIS is
%                    the density that enters buoyancy/Archimedes-type force
%                    terms (K_phi, K_eff) downstream — not phi, not rho_particle.
%     phi          : packing fraction = rho_bulk / rho_particle.
%
%   Conditions:
%     GB   -> full | shallow      (bed depth; always compacted packing)
%     CHIN -> as_poured | dense    (packing state; full container only)
%
%   If (material, condition) is unknown, returns NaNs with ok=false so callers
%   can flag rather than silently mis-scale.

    material  = upper(strtrim(char(material)));
    condition = lower(strtrim(char(condition)));

    % {material, condition, rho_particle_g_cm3, rho_bulk_g_cm3, phi}
    %   FINALISED VALUES. phi is the mean of five independent preparations
    %   (pour -> level -> weigh); rho_bulk = phi * rho_particle. These supersede
    %   the earlier single-measurement figures (0.629 / 0.636 / 0.276 / 0.409),
    %   which came from one preparation each and are no longer used anywhere.
    %   The +/- on phi (0.004 / 0.009 / 0.004 / 0.004) is RANDOM only; particle
    %   densities are treated as exact, so absolute phi carries an additional
    %   systematic that cancels between conditions but not against literature.
    T = {
        'GB',   'full',      2.50, 1.560, 0.624
        'GB',   'shallow',   2.50, 1.607, 0.643
        'CHIN', 'as_poured', 2.35, 0.658, 0.280
        'CHIN', 'dense',     2.35, 0.945, 0.402 };

    sub = struct('material',material, 'condition',condition, ...
                 'rho_particle_g_cm3',NaN, 'rho_bulk_g_cm3',NaN, 'phi',NaN, ...
                 'ok',false, 'reason','');

    idx = find(strcmpi(T(:,1), material) & strcmpi(T(:,2), condition), 1);
    if isempty(idx)
        sub.reason = sprintf('no substrate entry for (%s, %s)', material, condition);
        return;
    end
    sub.rho_particle_g_cm3 = T{idx,3};
    sub.rho_bulk_g_cm3     = T{idx,4};
    sub.phi                = T{idx,5};
    sub.ok                 = true;
end
