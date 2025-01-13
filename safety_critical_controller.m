function [u, scc_active] = safety_critical_controller(Pagc_in, Pload, Pdc, gamma , x_state, ...
    A, B, C)
    %{
        INPUTS: 
        - Pagc_in: float(4,1): is the agc control signal (pu) to each of areas 1, 2, 3, 4
        - Pload: float(4,1): is the load change (pu) in areas 1, 2, 3, 4
        - Pdc: float(4,1): is the dc power injection (pu) to areas 1, 2, 3, 4
        - x_state: float(12,1): is the state of the system [Pg (pu), frequency (pu), Ptie
        (pu)]1,2,3,4
        
        - alpha_scc: float(1,1): is the alpha parameter of the safety-critical controller.
        alpha_scc = 0 deactivates the safety critical controller.
        Lower alpha_scc makes the safety critical controller less relaxed.
        
        OUTPUTS:
        - u: float(12,1): is the input to the system consisting of [Pagc (pu), Pload (pu), 
        Pds (pu)]1,2,3,4
        - scc_active \in {0,1} is a boolean for whether the safety critical controller
        was activated or not.
        
    %}
    
    if gamma == 0
        scc_active = 0;
        u = [
            Pagc_in(1); Pload(1); Pdc(1);
            Pagc_in(2); Pload(2); Pdc(2);
            Pagc_in(3); Pload(3); Pdc(3);
            Pagc_in(4); Pload(4); Pdc(4)
        ];
    else
        options = optimoptions('fmincon', ...
            'Algorithm', 'sqp', 'MaxIterations', 500);
        
        Pagc_out = fmincon(...
            @(y) norm(Pagc_in-y), zeros(4,1), ...
            [], [], [], [], -0.3*ones(4,1), 0.3*ones(4,1), ...
            @(y) dhdt(y, Pload, Pdc, A, B, C, x_state, gamma), ...
            options ...
            );
        
        if Pagc_in ~= Pagc_out
            scc_active = 1;
        else
            scc_active = 0;
        end
        
        u = [
            Pagc_out(1); Pload(1); Pdc(1);
            Pagc_out(2); Pload(2); Pdc(2);
            Pagc_out(3); Pload(3); Pdc(3);
            Pagc_out(4); Pload(4); Pdc(4)
        ];
        
    end
end

function [c,ceq] = dhdt(u, PL, Pdc, A, B, C, x, gamma)
    %{
        u: float(4,1): is the Pagc
        x: float(12,1): is the state of the system
        alpha: float(1,1): is the SCC relaxation constant
        PL: float(4,1): is the load changes in the system
        Pdc: float(4,1): is the DC power injections
    %}
    Fsquared = (0.4/60)^2*ones(4,1); % 0.5 Hz to pu (4x1)

    frequency = C*x; % area frequencies in pu (4x1)
    
    h = max( eps, Fsquared - frequency.^2 ); % barrier function (4x1)
    
    Bl = -log( h ./ (1 + h) ); % logarithmic barrier function (4x1)

    dBldh = diag(-1 ./ (h + h.^2)); % 4x4

    dhdx = -2 * diag(frequency) * C; % 4x4 x 4x12 = 4x12

    LfB = dBldh * dhdx * A * x; % 4x4 x 4x12 x 12x12 x 12x1 = 4x1

    LgB = dBldh * dhdx * B(:, [1,4,7,10]); % 4x4 x 4x12 x 12x4 = 4x4

    LhB = dBldh * dhdx * (B(:, [2,5,8,11])*PL + B(:, [3,6,9,12])*Pdc);
    % 4x4 x 4x12 x 12x4 x 4x1 = 4x1
   
    c = LfB + LhB + LgB*u - gamma./Bl; % inequality constraint 4x1
    ceq = [];
end