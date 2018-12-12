function betas = plate_stress_regression(filename)
    
    % Import FEA results
    [num, raw] = xlsread(filename);
    thickness = [1.93;1.61;2.99;1.32;0.19;0.59;2.55;1.14;2.36;0.84];

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

    % Split data

    % shuffle the dataset
    sd = 1;
    rng(sd);
    rand_pos = randperm(height(raw)); % positions

    shuffled_raw= raw(rand_pos,:);

    splitPt1 = floor(0.75*height(raw));  % Case 1 split point

    training = shuffled_raw(1:splitPt1,:);
    test = shuffled_raw(splitPt1+1:end,:);

    % Regression Models using training data

    modelfun = @(b,x)(b(1)*exp(-b(2)*x + b(3))) + b(4);  % Model form
    beta0 = [2,6,8,30];  % Guesses
    opts = statset('nlinfit');
    opts.RobustWgtFun = 'bisquare';  % ignore outliers

    % Regression model for each material
    steelbeta = nlinfit(training.thickness, training.steel, modelfun, beta0);
    aluminiumbeta = nlinfit(training.thickness, training.aluminium, modelfun, beta0);
    zincbeta = nlinfit(training.thickness, training.zinc, modelfun, beta0);
    magnesiumbeta = nlinfit(training.thickness, training.magnesium, modelfun, beta0);
    
    % Return betas
    betas = [steelbeta; aluminiumbeta; zincbeta; magnesiumbeta];
end


