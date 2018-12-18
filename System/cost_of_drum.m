function [data] = cost_of_drum(betas)
    
    % Initiate table for results
    data = table();
    
    % Arrays for looping through materials/radii and their values
    materials = ["Aluminium","ABS","Copper","Zinc-Aluminium Alloy"];
    thicknesses = [0:0.1:1.5];
    price = [1.62 2.25 4.48 2.42];
    density = [2900 1213 8943 7000];
    yield_max = [27.5 35 34  35.9];
    modelfun = @(b,x)(b(1)*exp((-b(2)*x)+b(3))+b(4));

    % Looping through each combination of materials and radii 
    for x = 1:length(materials)
        for i = 1:length(thicknesses)
            
            % Finding the current values for each iteration
            current_material = materials(x);
            thickness_current = thicknesses(i);
            current_beta = betas(:,x);
            current_stress = modelfun(current_beta,thickness_current);

            % Calculating the volume for current iteration
            area = (pi()*(0.28+(thickness_current/1000))^2) - (pi()*0.28^2);
            cylinder = area*0.4;
            sides = (pi()*(0.28+(thickness_current/1000))^2)*(thickness_current/1000);
            volume = cylinder + sides;

            % Calculating the cost for current iteration
            price_current = price(x);
            density_current = density(x);
            mass = volume * density_current;
            cost = price_current * mass;

            % Calculating omega_d
            r_p = 0.0716;
            omega_p = 250;
            r_d = 0.28+(i/1000);
            omega_d = (r_p*omega_p)/r_d;

            % Add the results from current iteration to data array 
            new_row = {current_material, thickness_current, cost, omega_d};

            % But only the solutions that are successful
            if current_stress < yield_max(x)
                data = [data;new_row];

            end

        end
    end
end