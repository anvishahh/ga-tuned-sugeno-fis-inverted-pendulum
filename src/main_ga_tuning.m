%% GA-Tuned Sugeno FIS Controller for Inverted Pendulum (Portfolio Version)
% This script tunes a Sugeno FIS PD controller using a Genetic Algorithm (GA)
% and saves key plots to the /results folder.
%
% Requirements:
% - Fuzzy Logic Toolbox
% - Global Optimization Toolbox

clear; close all; clc; format compact

%% Reproducibility seed (replace with your own seed if you want)
rng(284631 + 276061 + 286715);

%% Constants (plant + controller settings)
constants = struct();
constants.M_CART = rand() * 0.1 + 2;    % cart mass, kg
constants.M_PEND = 1;                   % pend mass, kg
constants.DT = 0.01;                    % sim interval, s
constants.G = 9.8;
constants.CTRL_STEP = 5;                % controller period = DT*CTRL_STEP
constants.L_PEND = 1;                   % pend length, m
constants.KF_CART = 0.3;                % cart friction coeff
constants.KF_PEND = 0.4;                % pend friction coeff
constants.F_MAX = 50;                   % max force on cart, N
constants.PHI = (rand() * 2 + 30);      % slope angle, deg
constants.T_SET = 1;                    % settling time, s
constants.T_TOTAL = 5;                  % total sim time, s
constants.GAMMA_SET = 1;                % settling angle band, deg

constants.DT2 = constants.DT * constants.DT;
constants.S_PHI = sind(constants.PHI);
constants.C_PHI = cosd(constants.PHI);
constants.M_TOTAL = constants.M_CART + constants.M_PEND;

constants.L76 = -constants.L_PEND * 7.0 / 6.0;
constants.MML76 = constants.M_TOTAL * constants.L76;
constants.WEIGHT = constants.G * constants.S_PHI * constants.M_TOTAL;

constants.START_GAMMAS = [-10, -5, 5, 10]; % deg
constants.P_ERR_UNI = [-20, 20];
constants.D_ERR_UNI = [-40, 40];
constants.F_UNI = [-50, 50];

%% Build initial Sugeno FIS structure
fis = initFis(constants);

% Keep these in base workspace (objFcn uses evalin for simplicity)
assignin("base","fis",fis);
assignin("base","constants",constants);

%% GA setup
[ga_A, ga_b, ga_lb, ga_ub, ga_options] = initGAOptions(constants);
ga_options.OutputFcn = @outputFcn;

NUM_PARAMS = 25;
[~,~,ga_exit_flag,ga_output] = ga(@objFcn, NUM_PARAMS, ga_A, ga_b, [], [], ga_lb, ga_ub, [], ga_options); %#ok<ASGLU>

% Pull best history from base
history_best = evalin("base","history_best");
best_chrom = history_best{end,1};

%% Build tuned system and simulate
tuned_fis = adjustMFs(fis, best_chrom);

tuned_responses = [];
for k = 1:numel(constants.START_GAMMAS)
    tuned_responses = [tuned_responses, simFis(tuned_fis, constants, constants.START_GAMMAS(k))]; %#ok<AGROW>
end

%% FIG4: best fitness per generation
fig4 = figure(4);
set(fig4, 'Visible','on');
fig4.Position = [0, 0, 720, 240];
semilogy([history_best{:,2}]);
grid on
title("Fitness of Best Chromosome Across Generations");
ylabel("Fitness Score (lower is better)")
xlabel("Generation")

%% FIG2: membership functions + output surface
fig2 = figure(2);
set(fig2,'Visible','on')
fig2.Position = [0, 0, 1600, 300];
tld2 = tiledlayout(1,3,"TileSpacing","tight");
title(tld2,"Tuned Membership Functions + Output Surface")

nexttile
plotmf(tuned_fis,'input',1,1000)
title("p\_err MFs","Interpreter","none"); grid on

nexttile
plotmf(tuned_fis,'input',2,1000)
title("d\_err MFs","Interpreter","none"); grid on

nexttile
opts = gensurfOptions;
opts.NumGridPoints = 50;
gensurf(tuned_fis, opts);
title("Sugeno Output Surface","Interpreter","none"); grid on

%% FIG3: responses (gamma, gamma_vel, force)
fig3 = figure(3);
set(fig3,'Visible','on')
fig3.Position = [0, 0, 1200, 720];

tld3 = tiledlayout(numel(constants.START_GAMMAS), 3, "TileSpacing","tight");
title(tld3, "Tuned FIS Responses");
xlabel(tld3, "Time (s)");

for i = 1:numel(constants.START_GAMMAS)
    r = tuned_responses(i);
    res = evalResponse(constants, r);
    start_gamma = constants.START_GAMMAS(i);

    nexttile
    plot(r.time, r.gamma); grid on
    title(sprintf('\\gamma (start=%.1f deg)', start_gamma));
    xlim([-0.1, r.time(end)])
    if ~isnan(res.settling_time)
        xline(res.settling_time,'--k');
    end
    ylabel("\gamma (deg)")

    nexttile
    plot(r.time, r.gamma_vel); grid on
    title('\gamma velocity');
    xlim([-0.1, r.time(end)])
    ylabel("d\gamma/dt (deg/s)")

    nexttile
    plot(r.time, r.force); grid on
    title('Cart force');
    xlim([-0.1, r.time(end)])
    ylabel("F (N)")
end

%% Save figures
if ~exist("results","dir")
    mkdir("results");
end
exportgraphics(fig2,'results/MFs.png','BackgroundColor','white');
exportgraphics(fig3,'results/tuned.png','BackgroundColor','white');
exportgraphics(fig4,'results/history_best.png','BackgroundColor','white');

disp("Saved results to /results:");
disp(" - results/MFs.png");
disp(" - results/tuned.png");
disp(" - results/history_best.png");

disp("Tuned Sugeno output constants:");
disp(tuned_fis.Output.MembershipFunctions);

%% ========================= Local functions =========================

function fis = initFis(constants)
    fis = sugfis('Name', 'ga_tuned_pd_sugeno');

    fis = addInput(fis, constants.P_ERR_UNI, Name='p_err');
    fis = addInput(fis, constants.D_ERR_UNI, Name='d_err');
    fis = addOutput(fis, constants.F_UNI, Name='f');

    % 5 Gaussian MFs per input
    terms = ["NB","NS","ZO","PS","PB"];
    for t = 1:5
        fis = addMF(fis,'p_err','gaussmf',[1,0],Name=terms(t));
        fis = addMF(fis,'d_err','gaussmf',[1,0],Name=terms(t));
    end

    % 5 constant outputs (Sugeno)
    for t = 1:5
        fis = addMF(fis,'f','constant',0,Name=terms(t));
    end

    % PD fuzzy rule table (standard pattern)
    rule_table = [
        "NB","NB","NS","NS","ZO";
        "NB","NS","NS","ZO","PS";
        "NS","NS","ZO","PS","PS";
        "NS","ZO","PS","PS","PB";
        "ZO","PS","PS","PB","PB"
    ];

    rules = strings(0,1);
    for i = 1:5
        for j = 1:5
            rules(end+1,1) = sprintf("p_err==%s & d_err==%s => f==%s", ...
                terms(i), terms(j), rule_table(i,j)); %#ok<AGROW>
        end
    end
    fis = addRule(fis, rules);
end

function fis = adjustMFs(fis, chromosome)
    % Chromosome layout:
    % p_err: [std mean] x5 => 10 params
    % d_err: [std mean] x5 => 10 params
    % f constants: 5 params => total 25

    % p_err
    fis.Input(1).MembershipFunctions(1).Parameters = chromosome(1:2);
    fis.Input(1).MembershipFunctions(2).Parameters = chromosome(3:4);
    fis.Input(1).MembershipFunctions(3).Parameters = chromosome(5:6);
    fis.Input(1).MembershipFunctions(4).Parameters = chromosome(7:8);
    fis.Input(1).MembershipFunctions(5).Parameters = chromosome(9:10);

    % d_err
    fis.Input(2).MembershipFunctions(1).Parameters = chromosome(11:12);
    fis.Input(2).MembershipFunctions(2).Parameters = chromosome(13:14);
    fis.Input(2).MembershipFunctions(3).Parameters = chromosome(15:16);
    fis.Input(2).MembershipFunctions(4).Parameters = chromosome(17:18);
    fis.Input(2).MembershipFunctions(5).Parameters = chromosome(19:20);

    % Output constants
    fis.Output.MembershipFunctions(1).Parameters = chromosome(21);
    fis.Output.MembershipFunctions(2).Parameters = chromosome(22);
    fis.Output.MembershipFunctions(3).Parameters = chromosome(23);
    fis.Output.MembershipFunctions(4).Parameters = chromosome(24);
    fis.Output.MembershipFunctions(5).Parameters = chromosome(25);
end

function [p_err, d_err] = getErrors(gamma, gamma_vel)
    % PD errors (reference = 0)
    p_err = 0 - gamma;
    d_err = 0 - gamma_vel;
end

function response = simFis(fis, constants, start_gamma)
    T = 0:constants.DT:constants.T_TOTAL;
    CS = nan(size(T));
    CV = nan(size(T));
    PS = nan(size(T)); % gamma (deg later)
    PV = nan(size(T));
    F  = nan(size(T));

    f = 0;
    theta = deg2rad(start_gamma - constants.PHI);
    pv = 0; cs = 0; cv = 0;

    is_stable = true;

    for i = 1:numel(T)
        C_THETA = cos(theta);
        S_THETA = sin(theta);

        tmp_a = f - 0.5*pv*pv*S_THETA - constants.WEIGHT - cv*constants.KF_CART;
        tmp_b = -constants.G * (S_THETA*constants.C_PHI + C_THETA*constants.S_PHI) - pv*constants.KF_PEND;
        tmp_denom = constants.MML76 - 0.5*C_THETA*C_THETA;

        ca = (constants.L76*tmp_a + 0.5*C_THETA*tmp_b) / tmp_denom;
        pa = (-C_THETA*tmp_a + constants.M_TOTAL*tmp_b) / tmp_denom;

        % Integrate (Euler-ish)
        theta = theta + pv*constants.DT + 0.5*pa*constants.DT2;
        ps = theta + deg2rad(constants.PHI);
        pv = pv + pa*constants.DT;

        cs = cs + cv*constants.DT + 0.5*ca*constants.DT2;
        cv = cv + ca*constants.DT;

        % Controller update
        if mod(i-1, constants.CTRL_STEP) == 0
            [p_err, d_err] = getErrors(rad2deg(ps), rad2deg(pv));
            p_err = p_err + randn()*1e-4;
            d_err = d_err + randn()*1e-4;

            if p_err < constants.P_ERR_UNI(1) || p_err > constants.P_ERR_UNI(2) || ...
               d_err < constants.D_ERR_UNI(1) || d_err > constants.D_ERR_UNI(2)
                is_stable = false;
                break;
            end

            f = evalfis(fis, [p_err, d_err]);
            f = clip(f, -constants.F_MAX, constants.F_MAX);
        end

        CS(i) = cs; CV(i) = cv; PS(i) = ps; PV(i) = pv; F(i) = f;
    end

    PS = rad2deg(PS);
    PV = rad2deg(PV);

    response = struct();
    response.time = T;
    response.gamma = PS;
    response.gamma_vel = PV;
    response.force = F;
    response.cart = CS;
    response.cart_vel = CV;
    response.is_stable = is_stable;
end

function results = evalResponse(constants, response)
    if response.is_stable
        results.settling_time = response.time(end);
        for i = numel(response.time):-1:1
            g = response.gamma(i);
            if (g < -constants.GAMMA_SET || g > constants.GAMMA_SET)
                if i < numel(response.time)
                    results.settling_time = response.time(i+1);
                end
                break
            end
        end

        results.last_gamma = response.gamma(end);

        sgn = sign(response.gamma(1));
        results.peak_fraction_gamma = abs(min(response.gamma * sgn) / response.gamma(1));

        idx = find(response.time > 1, 1, 'first');
        results.oscillation_gamma = std(response.gamma(idx:end));

    else
        results.last_gamma = constants.P_ERR_UNI(2);
        results.settling_time = constants.T_TOTAL;
        results.peak_fraction_gamma = abs(constants.P_ERR_UNI(2) / response.gamma(1));
        results.oscillation_gamma = constants.P_ERR_UNI(2);
    end
end

function total_score = objFcn(chromosome)
    fis       = evalin("base","fis");
    constants = evalin("base","constants");

    fis = adjustMFs(fis, chromosome);

    % Weights
    k_s   = 1.0;
    k_f   = 2.0;
    k_osc = 1.2;
    k_u   = 0.6;
    k_du  = 1.0;

    % Normalizers
    Tnorm = max(constants.T_TOTAL, 1e-9);
    Gnorm = max(constants.GAMMA_SET, 1e-9);
    Fnorm = max(constants.F_MAX, 1e-9);

    total_score = 0;

    for start_gamma = constants.START_GAMMAS
        r   = simFis(fis, constants, start_gamma);
        res = evalResponse(constants, r);

        n_ts  = res.settling_time / Tnorm;
        n_gf  = abs(res.last_gamma) / Gnorm;
        n_osc = res.oscillation_gamma / Gnorm;

        jj = (r.time > 1) & isfinite(r.force);
        Fu = r.force(jj);

        if ~isempty(Fu)
            rmsF = sqrt(mean(Fu.^2));
            n_u  = rmsF / Fnorm;

            dF    = diff(Fu);
            rmsdF = sqrt(mean(dF.^2));
            n_du  = rmsdF / Fnorm;
        else
            n_u = 1; n_du = 1;
        end

        score = k_s*n_ts^2 + k_f*n_gf^2 + k_osc*n_osc^2 + k_u*n_u^2 + k_du*n_du^2;
        total_score = total_score + score;
    end
end

function [ga_A, ga_b, ga_lb, ga_ub, ga_options] = initGAOptions(constants)
    NUM_PARAMS = 25;
    NUM_POP = 100;
    THRESH = 1e-3;

    ga_lb = zeros(NUM_PARAMS,1);
    ga_lb(2:2:10)  = constants.P_ERR_UNI(1);
    ga_lb(12:2:20) = constants.D_ERR_UNI(1);
    ga_lb(1:2:19)  = 0; % std >= 0
    ga_lb(21:25)   = constants.F_UNI(1);
    ga_lb = ga_lb + THRESH;

    ga_ub = zeros(NUM_PARAMS,1);
    ga_ub(2:2:10)  = constants.P_ERR_UNI(2);
    ga_ub(12:2:20) = constants.D_ERR_UNI(2);
    ga_ub(1:2:19)  = Inf; % std unbounded
    ga_ub(21:25)   = constants.F_UNI(2);
    ga_ub = ga_ub - THRESH;

    % Linear inequality constraints to enforce mean ordering:
    % mean_NB < mean_NS < mean_ZO < mean_PS < mean_PB (for both inputs)
    ga_A = [
        0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0;

        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1;
    ];
    ga_b = zeros(size(ga_A,1),1);

    ga_options = optimoptions("ga");
    ga_options.PopulationSize = NUM_POP;
    ga_options.MaxGenerations = 50;
    ga_options.PlotFcn = {'gaplotscorediversity','gaplotbestf'};

    % Initial population
    population = zeros(NUM_POP, NUM_PARAMS);
    for i = 1:NUM_POP
        means_p = sort(unifrnd(constants.P_ERR_UNI(1), constants.P_ERR_UNI(2), 1, 5));
        means_d = sort(unifrnd(constants.D_ERR_UNI(1), constants.D_ERR_UNI(2), 1, 5));
        f = sort(unifrnd(constants.F_UNI(1)+THRESH, constants.F_UNI(2)-THRESH, 1, 5));

        params_p = [1, means_p(1), 1, means_p(2), 1, means_p(3), 1, means_p(4), 1, means_p(5)];
        params_d = [1, means_d(1), 1, means_d(2), 1, means_d(3), 1, means_d(4), 1, means_d(5)];

        population(i,:) = [params_p, params_d, f];
    end
    ga_options.InitialPopulationMatrix = population;
end

function [state,options,optchanged] = outputFcn(options,state,flag)
    persistent history_best
    optchanged = false;

    ibest = state.Best(end);
    ibest = find(state.Score == ibest, 1, 'last');
    bestx = state.Population(ibest,:);
    bestf = objFcn(bestx);

    switch flag
        case 'init'
            history_best = {bestx, bestf};
        case {'iter','done'}
            if state.Generation > size(history_best,1)
                history_best = [history_best; {bestx, bestf}]; %#ok<AGROW>
            end
    end

    assignin('base','history_best',history_best);
    fprintf("\nGen %d Best Obj=%f\n", state.Generation, bestf);
end

function y = clip(x, lo, hi)
    y = min(max(x, lo), hi);
end
