tic
clc
close all
clear all

% Find cost combinations of plate
filename = 'subsystem1_results.xlsx';
betas = plate_stress_regression(filename);
all_plate_options = cost_of_plate(betas);

% Find cost of combinations drum
filename = 'subsystem2_results.csv';
betas = drum_stress_regression(filename);
all_drum_options = cost_of_drum(betas);

% Maximum spend is �1 
prices_plate = all_plate_options{:,3};
prices_drum = all_drum_options{:,3};

% Find all options under �3
total_costs = [];
for p = 1:length(prices_plate)
    for d = 1:length(prices_drum)
        total_cost = prices_drum(d) + prices_plate(p);
        
        drum = all_drum_options{d,:};
        plate = all_plate_options{p,:};

        new_row = [drum plate total_cost];
        
        if total_cost < 3 % Must be under �3
            total_costs = [total_costs; new_row];
        end
       
    end
end

% normalise drum speed and heat flux
% [red_wine_x_newTrain1,PS_rwxTrain1] = mapstd(red_wine_x_new(1:splitPt1,:)');
heat_flux = str2double(total_costs(:,8));
drum_v = str2double(total_costs(:,4));
heat_flux_scaled = rescale(heat_flux);
drum_v_scaled = rescale(drum_v);

weights = [];
% Find option with best drying efficiency
for r = 1:length(total_costs)
    heat_flux_current = heat_flux_scaled(r);
    drum_speed_current = drum_v_scaled(r);
    
    %Weights
    weight = -0.5*heat_flux_current + 0.5*drum_speed_current;
    weights = [weights; weight];
    %store

end

all = [total_costs weights];

sorted = sortrows(all, 10);
top = sorted(1,:);
headings = ["Drum Material", "Drum thickness (mm)", "Drum Cost (�)", "Drum Speed (RPM)", "Plate Material", "Plate Thickness (m)", "Plate Cost (�)", "Heat Flux (W/m2)", "Total Cost (�)", "Weighting"];
optimum = [headings; top]
toc