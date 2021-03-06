function [multi_objective_pass] = cost_of_plate(betas)
    % Function to return all possible plate costs that pass FoS test

    % Data to be used
    % All in order steel - aluminium - zinc - magnesium

    % Yield Strengths
    yield_strengths = [351.517, 55.1458, 220, 148];  %MPa

    % Thermal conductivities
    ks = [47, 170, 108.9, 160]; %W/mK

    % material cost per unit (average upper and lower bounds)
     % GBP per kg
    steel_cost = (0.531 + 0.577)/2; % GBP per kg
    aluminium_cost = (1.51 + 1.62)/2;
    zinc_cost = (1.75 + 2.42)/2;
    magnesium_cost = (1.9 + 2.12)/2;
    costs = [steel_cost, aluminium_cost, zinc_cost, magnesium_cost];

    % material densities
    densities = [7900, 2700, 6700, 1700]; %kg/m3

    % Empty arrays
    multi_objective_pass = table();
    multi_objective_fail = table();

    % Optimised parameters from subsystem optimisation
    Theater = 680.15;
    v = 15;
    LAir = 0;
    r = 0.225;  % m
    thickness_range = 0:0.1:3;
    safety_factor = 2;

    materials = ["Steel", "Aluminium", "Zinc", "Magnesium"];
    modelfun = @(b,x)(b(1)*exp(-b(2)*x + b(3))) + b(4);  % Model form
    
    
    for x = 1:length(thickness_range)
        for i = 1:length(materials)
            % Find stress
            current_beta = betas(i,:);
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
            
            % Check yield strength
            if  current_stress < safety_factor * yield_strengths(i)
                multi_objective_pass = [multi_objective_pass; newrow];

            else
                multi_objective_fail = [multi_objective_fail; newrow];

            end  
        end
    end

    
end

%% Finding heatflux in multi-objective optimisation
function q = heatfluxpairs(x, k, Theater, v, LAir) % x in mm
        
            % Parameters
            Tinside = 330.4;
            kair = 2.816;
           
            q = -((Theater-Tinside)/(1/h(v) + LAir/kair + 1/h(v) + (x/1000)/k + 1/h(v)));
        
end

%% Heat transfer coeffcient
function h_c = h(v)
    h_c = 10.45-v+10*v^0.5;
end