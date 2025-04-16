function [Best_F, Best_P, Conv_curve] = COVIDOA(nPop, MaxIt, minVal, maxVal, D, CostFunction, MR, shifttingNo, numOfSubprotiens)
    VarMin = ones(1, D) .* minVal;
    VarMax = ones(1, D) .* maxVal;
    gamma = 0.5;
    beta = 0.5;
    bestsol.Cost = inf;
    bestsol.Position = zeros(1, D);
    empty_individual.Position = [];
    empty_individual.Cost = [];
    pop = repmat(empty_individual, nPop, 1);

    % Initialization
    for i = 1:nPop
        pop(i).Position = unifrnd(VarMin, VarMax, [1 D]);
        pop(i).Cost = CostFunction(pop(i).Position);
        if pop(i).Cost < bestsol.Cost
            bestsol.Cost = pop(i).Cost;
            bestsol.Position = pop(i).Position;
        end
    end

    Conv_curve = nan(MaxIt, 1);
    noImprovementCount = 0;
    restartThreshold = 20; % Number of iterations with no improvement before restart
    originalMR = MR;

    % Main Loop
    for it = 1:MaxIt
        % Adaptive parameters
        gamma = 0.5 + 0.5 * (1 - it / MaxIt); % Example of dynamic adjustment
        beta = 0.5 + 0.5 * (it / MaxIt); % Example of dynamic adjustment

        % Virus replication phase through frameshifting technique
        c = [pop.Cost];
        avgc = mean(c);
        if avgc ~= 0
            c = c / avgc;
        end
        probs = exp(-beta * c);
        x = [];
        for k = 1:nPop
            parent = pop(TournamentSelection(probs, 3)); % Using tournament selection
            parent.Position = max(parent.Position, VarMin);
            parent.Position = min(parent.Position, VarMax);
            x(k,:) = parent.Position;
            for t = 1:numOfSubprotiens
                for i = 1:D - shifttingNo
                    x(k,i) = x(k,i + shifttingNo);
                end
                r = unifrnd(minVal, maxVal, [1 shifttingNo]);
                x(k,:) = [x(k,1:D - shifttingNo) r];
                subprotien(t,:) = x(k,:);
            end
            newvirus(k,:) = UniformCrossover(subprotien(1,:), subprotien(2,:), gamma);
            newvirus(k,:) = max(newvirus(k,:), VarMin);
            newvirus(k,:) = min(newvirus(k,:), VarMax);
        end

        for t = 1:nPop
            childcost = CostFunction(newvirus(t,:));
            if childcost < bestsol.Cost
                bestsol.Position = newvirus(t,:);
                bestsol.Cost = childcost;
                noImprovementCount = 0; % Reset counter if improvement
            end
        end

        newPop = repmat(empty_individual, nPop, 1);
        for l = 1:nPop
            for k = 1:D
                R = rand();
                if R < MR
                    newvirus(l,k) = minVal + rand * (maxVal - minVal);
                end
                newPop(l).Position = newvirus(l,:);
            end
            newPop(l).Position = max(newPop(l).Position, VarMin);
            newPop(l).Position = min(newPop(l).Position, VarMax);
            newPop(l).Cost = CostFunction(newPop(l).Position);
            if newPop(l).Cost < bestsol.Cost
                bestsol.Position = newPop(l).Position;
                bestsol.Cost = newPop(l).Cost;
                noImprovementCount = 0; % Reset counter if improvement
            end
        end

        % Implement elitism
        pop = SortPopulation([pop; newPop]);
        pop = pop(1:nPop);

        Conv_curve(it) = bestsol.Cost;
        Best_F = bestsol.Cost;
        Best_P = bestsol.Position;

        disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(Conv_curve(it))]);

        % Check for no improvement
        if it > 1 && Conv_curve(it) >= Conv_curve(it - 1)
            noImprovementCount = noImprovementCount + 1;
        else
            noImprovementCount = 0;
        end

        % Adaptive Mutation Rate
        if noImprovementCount >= 10
            MR = min(MR * 1.5, 0.5); %Increase mutation rate, cap at 0.5
        else
            MR = originalMR; % Reset to original mutation rate
        end

        % Population Restart
        if noImprovementCount >= restartThreshold
            disp('Restarting population due to stagnation...');
            for i = 1:nPop
                pop(i).Position = unifrnd(VarMin, VarMax, [1 D]);
                pop(i).Cost = CostFunction(pop(i).Position);
                if pop(i).Cost < bestsol.Cost
                    bestsol.Cost = pop(i).Cost;
                    bestsol.Position = pop(i).Position;
                end
            end
            noImprovementCount = 0; % Reset counter after restart
        end

        % Hybridize with a local search
        if mod(it, 10) == 0
            bestsol = LocalSearch(bestsol, CostFunction, VarMin, VarMax);
        end
    end
end

% Tournament Selection Function
function i = TournamentSelection(probs, t_size)
    % Randomly select t_size individuals and return the one with the best (lowest) cost
    selected = randi(length(probs), [1, t_size]);
    [~, idx] = min(probs(selected));
    i = selected(idx);
end

% Local Search Function
function bestsol = LocalSearch(bestsol, CostFunction, VarMin, VarMax)
    % Perform a simple local search around the best solution
    step = (VarMax - VarMin) * 0.01;
    for i = 1:length(bestsol.Position)
        newsol = bestsol;
        newsol.Position(i) = newsol.Position(i) + step(i);
        newsol.Position = max(newsol.Position, VarMin);
        newsol.Position = min(newsol.Position, VarMax);
        newsol.Cost = CostFunction(newsol.Position);
        if newsol.Cost < bestsol.Cost
            bestsol = newsol;
        end
    end
end
