function [data] = cost_of_drum(betas)
    data = table();

    materials = ["Aluminium","ABS","Copper","Zinc-Aluminium Alloy"];
    radii = [0:0.1:1.5];

    price = [1.62 2.25 4.48 2.42];
    density = [2900 1213 8943 7000];
    yield_max = [27.5 35 34  35.9];
    modelfun = @(b,x)(b(1)*exp((-b(2)*x)+b(3))+b(4));

    for x = 1:length(materials)
        for i = 1:length(radii)
            current_material = materials(x);
            radius_current = radii(i);
            current_beta = betas(:,x);
            current_stress = modelfun(current_beta,radius_current);

            area = (pi()*(0.28+(radius_current/1000))^2) - (pi()*0.28^2);
            cylinder = area*0.4;
            sides = (pi()*(0.28+(radius_current/1000))^2)*(radius_current/1000);
            volume = cylinder + sides;

            price_current = price(x);
            density_current = density(x);
            mass = volume * density_current;
            cost = price_current * mass;

            r_p = 0.0716;
            omega_p = 250;
            r_d = 0.28+(i/1000);
            omega_d = (r_p*omega_p)/r_d;

            % Add data to 
            new_row = {current_material, radius_current, cost, omega_d};


            if current_stress < yield_max(x)
                data = [data;new_row];

            end

        end
    end
end