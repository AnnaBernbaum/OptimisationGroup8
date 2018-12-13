clc
close all

%% Data to be used
% Yield Strengths
yield_strength_zinc = 220; % MPa
yield_strength_steel = 351.571; % Mpa
yield_strength_aluminium = 55.1458; % MPa
yield_strength_magnesium = 148; % MPa
yield_strengths = [yield_strength_steel, yield_strength_aluminium, yield_strength_zinc, yield_strength_magnesium];

% Thermal conductivities
k_zinc = 108.9;  %W/mK
k_steel = 47;
k_aluminium = 170;
k_magnesium = 160;
ks = [k_steel, k_aluminium, k_zinc, k_magnesium];

% material cost per unit (average upper and lowed bounds)
 % GBP per kg
steel_cost = (0.531 + 0.577)/2; % GBP per kg
aluminium_cost = (1.51 + 1.62)/2;
zinc_cost = (1.75 + 2.42)/2;
magnesium_cost = (1.9 + 2.12)/2;
costs = [steel_cost, aluminium_cost, zinc_cost, magnesium_cost];

% material densities
steel_density = 7900; % kg/m3
aluminium_density = 2700;
zinc_density = 6700;
magnesium_density = 1700;
densities = [steel_density, aluminium_density, zinc_density, magnesium_density];

materials = ["Steel", "Aluminium", "Zinc", "Magnesium"];

global Tinside % Kelvin = 30 degrees celsius;
Tinside = 330.4;
global kair;% W/m K
kair = 2.816;

%% Sub-subsystem Optimisation

% Sampling: Latin hypercube
datapoints = lhsdesign(10,1, 'criterion', 'maximin');  % 10 samples for each material
datapoints = datapoints*3;  % scale to be between 0 and 3 mm

% as the results vary every time, the points used were: 
thickness = [1.93;1.61;2.99;1.32;0.19;0.59;2.55;1.14;2.36;0.84];

% Import FEA results
filename = 'Results.xlsx';
[num, raw] = xlsread(filename);

% Empty arrays for data
steel = zeros(10,1);
aluminium = zeros(10,1);
zinc = zeros(10,1);
magnesium = zeros(10,1);

j = 1;  % iterative to move through thicknesses

% read values from generated CSV
for i = 1:length(raw)
    
    stress = string(raw(3, i));
    
    stress = sscanf(stress, '%f');  % remove units, convert to double
    
    if string(raw(2, i)) == 'AISI 1020 @SOLIDWORKS Materials'
        steel(j) = stress;
    end
    if string(raw(2, i)) == '6061 Alloy @SOLIDWORKS Materials'
        aluminium(j) = stress;
    end
    if string(raw(2, i)) == 'Zinc AC41A Alloy, As Cast @SOLIDWORKS Materials'
        zinc(j) = stress;
    end
    if string(raw(2, i)) == 'Magnesium Alloy @SOLIDWORKS Materials'
        magnesium(j) = stress;
    end
    
    j = j+1;
    
    if j > 10  % Reset J after the 10 results
        j = 1;
    end
     
end

% Create table of results
raw = table(thickness, steel, aluminium, zinc, magnesium);

%% Split data

% shuffle the dataset
sd = 1;
rng(sd);
rand_pos = randperm(height(raw)); % positions

shuffled_raw= raw(rand_pos,:);

splitPt1 = floor(0.75*height(raw));  % Case 1 split point

training = shuffled_raw(1:splitPt1,:);
test = shuffled_raw(splitPt1+1:end,:);

%% Regression Models using training data

modelfun = @(b,x)(b(1)*exp(-b(2)*x + b(3))) + b(4);  % Model form
beta0 = [2,6,8,30];  % Guesses
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';  % ignore outliers

% Regression model for each material
steelbeta = fitnlm(training.thickness, training.steel, modelfun, beta0)
aluminiumbeta = fitnlm(training.thickness, training.aluminium, modelfun, beta0)
zincbeta = fitnlm(training.thickness, training.zinc, modelfun, beta0)
magnesiumbeta = fitnlm(training.thickness, training.magnesium, modelfun, beta0)

% Find optimal thickness of each metal
safety_factor = 2;

% rearranged equation such that x = 
reversemodelfun = @(b,f)((log((f - b(4))/b(1))-b(3))/(-b(2)));

% Find minimum value of x to satisfy constraint g5
t_steel = reversemodelfun(steelbeta.Coefficients.Estimate, (yield_strength_steel*2));
t_aluminium = reversemodelfun(aluminiumbeta.Coefficients.Estimate, yield_strength_aluminium*2);
t_zinc = reversemodelfun(zincbeta.Coefficients.Estimate, yield_strength_zinc*2);
t_magnesium = reversemodelfun(magnesiumbeta.Coefficients.Estimate, yield_strength_magnesium*2);
thicknesses = [t_steel, t_aluminium, t_zinc, t_magnesium];

% Calculate ratios
ratio_steel = (t_steel/1000)/k_steel; % minimise this for each metal
ratio_aluminium = (t_aluminium/1000)/k_aluminium;
ratio_zinc = (t_zinc/1000)/k_zinc; 
ratio_magnesium = (t_magnesium/1000)/k_magnesium; 

% Plot regression models and training and test data
x = 0:0.01:3;
figure(1)

subplot(2,2,1)
scatter(training.thickness, training.steel)  % Training points
hold on 
scatter(test.thickness, test.steel)  % Test points
hold on
plot(x, modelfun(steelbeta.Coefficients.Estimate, x))  % Model
hold on
plot(linspace(0,3), ones(100)*2*yield_strength_steel, 'g')
title('Steel')
xlabel('Thickness (mm)')
ylabel('Maximum von Mises Stress (MPa)')
legend('Training', 'Test', 'Model', 'Steel Yield Strength')
axis([0 3 0 1200])

subplot(2,2,2)
scatter(training.thickness, training.aluminium)
hold on 
scatter(test.thickness, test.aluminium)
hold on
plot(x, modelfun(aluminiumbeta.Coefficients.Estimate, x))
hold on
plot(linspace(0,3), ones(100)*2*yield_strength_aluminium, 'g')
title('Aluminium')
xlabel('Thickness (mm)')
ylabel('Maximum von Mises Stress (MPa)')
legend('Training', 'Test', 'Model', 'Aluminium Yield Strength')

subplot(2,2,3)
scatter(training.thickness, training.zinc)
hold on 
scatter(test.thickness, test.zinc)
hold on
plot(x, modelfun(zincbeta.Coefficients.Estimate, x))
hold on
plot(linspace(0,3), ones(100)*2*yield_strength_zinc, 'g')
title('Zinc')
xlabel('Thickness (mm)')
ylabel('Maximum von Mises Stress (MPa)')
legend('Training', 'Test', 'Model', 'Zinc Yield Strength')

subplot(2,2,4)
scatter(training.thickness, training.magnesium)
hold on 
scatter(test.thickness, test.magnesium)
hold on
plot(x, modelfun(magnesiumbeta.Coefficients.Estimate, x))
hold on
plot(linspace(0,3), ones(100)*2*yield_strength_magnesium, 'g')
title('Magnesium')
xlabel('Thickness (mm)')
ylabel('Maximum von Mises Stress (MPa)')
legend('Training', 'Test', 'Model', 'Magnesium Yield Strength')


%% Subsystem Optimisation

% set initial guesses
Theaterguess = 600; % Degrees K
vguess = 10; % m/s
LAirguess = 0.04; %m
x0 = [Theaterguess vguess LAirguess];

% Bounds
lb = [-Inf, 0, 0];
ub = [680.15, 15, Inf];

%%%% Fmincon optimisation with interior point %%%%
options = optimoptions('fmincon', 'Algorithm','interior-point');
tic  % start timer
[x_opt_steel,q_steel] = optimise_fmincon(ratio_steel, x0, lb, ub, options);
[x_opt_aluminium,q_aluminium] = optimise_fmincon(ratio_aluminium, x0, lb, ub, options);
[x_opt_zinc,q_zinc] = optimise_fmincon(ratio_zinc, x0, lb, ub, options);
[x_opt_magnesium,q_magnesium] = optimise_fmincon(ratio_magnesium, x0, lb, ub, options);
time_fmincon_ip = toc;  % end timer

% Pick minimum q from the 4 materials
[q_opt_fmincon_ip, I_ip] = min([q_steel, q_aluminium, q_zinc, q_magnesium]);
opt_material_fmincon_ip = materials(I_ip); % Optimum material
x_opts = [x_opt_steel,; x_opt_aluminium; x_opt_zinc; x_opt_magnesium]; 
x_opt_ip = x_opts(I_ip,:); % Optimum Theater, v and Lair

%%%% Fmincon optimisation with SQP %%%%
options = optimoptions('fmincon', 'Algorithm','sqp');
tic
[x_opt_steel,q_steel] = optimise_fmincon(ratio_steel, x0, lb, ub, options);
[x_opt_aluminium,q_aluminium] = optimise_fmincon(ratio_aluminium, x0, lb, ub, options);
[x_opt_zinc,q_zinc] = optimise_fmincon(ratio_zinc, x0, lb, ub, options);
[x_opt_magnesium,q_magnesium] = optimise_fmincon(ratio_magnesium, x0, lb, ub, options);
time_fmincon_sqp = toc;

% Pick minimum q from the 4 materials
[q_opt_fmincon_sqp, I_sqp] = min([q_steel, q_aluminium, q_zinc, q_magnesium]);
opt_material_fmincon_sqp = materials(I_sqp); % Store answer
x_opts = [x_opt_steel,; x_opt_aluminium; x_opt_zinc; x_opt_magnesium];
x_opt_sqp = x_opts(I_sqp,:);

%%%% GA optimisation %%%%

% Bounds
lb = [-Inf, 0, 0, 1];
ub = [680.15, 15, Inf, 4];
ratios = [ratio_steel, ratio_aluminium, ratio_zinc, ratio_magnesium];

tic % start timer
[x_opt_ga,q_opt_ga] = optimise_ga(ratios, lb, ub) ;
time_ga = toc; % end timer

opt_material_ga = materials(x_opt_ga(4)); % Find optimum material

%%%% Check it is a Global optimum %%%%
lb = [-Inf, 0, 0];
ub = [680.15, 15, Inf];

options = 'sqp';
[x_opt_steel_global,q_steel_global] = optimise_fmincon_global(ratio_steel, x0, lb, ub, options);
[x_opt_aluminium_global,q_aluminium_global] = optimise_fmincon_global(ratio_aluminium, x0, lb, ub, options);
[x_opt_zinc_global,q_zinc_global] = optimise_fmincon_global(ratio_zinc, x0, lb, ub, options);
[x_opt_magnesium_global,q_magnesium_global] = optimise_fmincon_global(ratio_magnesium, x0, lb, ub, options);


%% Multiobjective for cost v material

% Empty arrays
multi_objective_pass = table();
multi_objective_fail = table();

% Parameters
Theater = x_opt_sqp(1);
v = x_opt_sqp(2);
LAir = x_opt_sqp(3);
r = 0.225;  % m
thickness_range = 0:0.1:1;
betas = [steelbeta.Coefficients.Estimate, aluminiumbeta.Coefficients.Estimate, zincbeta.Coefficients.Estimate, magnesiumbeta.Coefficients.Estimate];

for x = 1:length(thickness_range)
    for i = 1:length(materials)
        
        % Find stress
        current_beta = betas(:,i);
        current_stress = modelfun(current_beta, thickness_range(x));
        
        % Find heat flux    
        heatflux = heatfluxpairs(thickness_range(x), ks(i), Theater, v, LAir);
        
        % Find Cost
        current_material = materials(i);  
        current_thickness = thickness_range(x)/1000; %m
        material_volume = current_thickness * pi * r^2;
        material_density = densities(i);
        material_mass = material_density * material_volume;
        material_cost = material_mass * costs(i);
        
        % Append data
        newrow = {current_material, current_thickness, material_cost, heatflux};

        if  current_stress < safety_factor * yield_strengths(i)
            multi_objective_pass = [multi_objective_pass; newrow];

        else
            multi_objective_fail = [multi_objective_fail; newrow];
            
        end  
    end
end

% Plot multiobjective problem
figure(2)
scatter(multi_objective_pass{:,4}, multi_objective_pass{:,3})  % Cost v heatflux
hold on
scatter(multi_objective_fail{:,4}, multi_objective_fail{:,3})
hold on
plot(-3984.657,0.2717,'k.','MarkerSize',20);  % Highlight Pareto Set
grid on
ylabel('Cost (£)')
xlabel('Heat Flux (W/m^2)')
title('Multi-Objective Optimisation')
legend('Pass', 'Fail')


%%%%% Find % improvement %%%%%%%
thickness_original = 0.001;
ratio = thickness/k_steel;

q_original = -((Theaterguess-Tinside)/(1/h(v) + LAirguess/kair + 1/h(vguess) + ratio + 1/h(vguess)));
fmincon_improvement_ip = ((abs(q_opt_fmincon_ip)-abs(q_original))/abs(q_original))*100;
fmincon_improvement_sqp = ((abs(q_opt_fmincon_sqp)-abs(q_original))/abs(q_original))*100;
ga_improvement = ((abs(q_opt_ga)- abs(q_original))/abs(q_original))*100;

original_vol = thickness_original * pi * r^2;
original_mass = steel_density * original_vol;
cost_original = original_mass * steel_cost;
cost_improvement = ((cost_original)-((0.2717))/(cost_original))*100;


% Show table of results
Method = ["fmincon (interior point)" ;"fmincon (SQP)"; "GA"];
opt_q = [q_opt_fmincon_ip; q_opt_fmincon_sqp; q_opt_ga];
opt_material = [opt_material_fmincon_ip; opt_material_fmincon_sqp; opt_material_ga];
opt_thickness = [thicknesses(I_ip); thicknesses(I_sqp);thicknesses(x_opt_ga(4))];
opt_velocity = [x_opt_ip(2); x_opt_sqp(2); x_opt_ga(2)];
opt_lair = [x_opt_ip(3); x_opt_sqp(3); x_opt_ga(3)];
opt_theater = [x_opt_ip(1); x_opt_sqp(1); x_opt_ga(1)];

comp_time = [time_fmincon_ip;time_fmincon_sqp; time_ga];
improvement = [fmincon_improvement_ip;fmincon_improvement_sqp; ga_improvement];

table(Method, opt_q, opt_material, opt_thickness, comp_time, improvement, opt_lair, opt_theater, opt_velocity)


%% Heat transfer coeffcient
function h_c = h(v)
    h_c = 10.45-v+10*v^0.5;
end


%% Fmincon algorithm optimisation
function [x,heat_flux] =  optimise_fmincon(ratio, x0, lb, ub, options) 
    % Return optimum x (vector of 3 variables Theater, v, Lair)
    [x,heat_flux] = fmincon(@heatflux, x0, [], [], [], [], lb, ub ,[], options);
    
    % Nested function that computes the objective function     
        function q = heatflux(x)
            % Parameters
            global Tinside
            global kair
            Theater = x(1);
            v = x(2);
            LAir = x(3);

            q = -((Theater-Tinside)/(1/h(v) + LAir/kair + 1/h(v) + ratio + 1/h(v)));
        
        end
end

%% Global fmincon optimisation
function [x,heat_flux] =  optimise_fmincon_global(ratio, x0, lb, ub, options) 
    % Return optimum x (vector of 3 variables Theater, v, Lair)
    opts = optimoptions(@fmincon,'Algorithm', options);
    problem = createOptimProblem('fmincon','x0',x0,'objective', @heatflux,'lb',lb,'ub',ub, 'options',opts);
    [x, heat_flux] = run(GlobalSearch,problem);
    
    % Nested function that computes the objective function     
        function q = heatflux(x)
            % Parameters
            global Tinside
            global kair
            Theater = x(1);
            v = x(2);
            LAir = x(3);

            q = -((Theater-Tinside)/(1/h(v) + LAir/kair + 1/h(v) + ratio + 1/h(v)));
        
        end
end

%% Genetic algorithm optimisation
function [x,heat_flux] =  optimise_ga(ratios, lb, ub) 
    % Return optimum x (vector of 3 variables Theater, v, Lair)
    [x,heat_flux] = ga(@heatflux, 4, [], [], [], [], lb, ub ,[], 4); 

    
    % Nested function that computes the objective function     
        function q = heatflux(x)
        
            % Parameters
            global Tinside
            global kair
            Theater = x(1);
            v = x(2);
            LAir = x(3);
            
            x(4) = ratios(x(4)); % map integer values to 
            ratio = x(4);

            q = -((Theater-Tinside)/(1/h(v) + LAir/kair + 1/h(v) + ratio + 1/h(v)));
        
        end
end

%% Finding heatflux in multi-objective optimisation
function q = heatfluxpairs(x, k, Theater, v, LAir) % x in mm
        
            % Parameters
            global Tinside
            global kair
           
            q = -((Theater-Tinside)/(1/h(v) + LAir/kair + 1/h(v) + (x/1000)/k + 1/h(v)));
        
end