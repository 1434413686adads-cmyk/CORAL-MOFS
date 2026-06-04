classdef CORAL < ALGORITHM
% <multi> <binary> <large/none> <constrained/none>
% CORAL: A Rank-Memory Search Framework for Multi-Objective Feature Selection
%
% Parameters:
%   B     --- 5   --- number of sparsity layers in stratified environmental selection
%   Nm    --- 3   --- number of elite individuals used for local refinement
%   rho   --- 0.4 --- rank-memory historical/current fusion coefficient
%   gamma --- 0.4 --- rank-memory bias strength

    properties(Access = private)
        freq
    end

    methods
        function main(Algorithm, Problem)
            [B, Nm, rho, gamma] = Algorithm.ParameterSet(5, 3, 0.4, 0.4);

            F        = 0.5;
            CR       = 0.9;
            sigma    = 0.1;
            Fk       = 0.3;
            Kmin     = 1;
            KmaxFrac = 0.2;

            memN = Nm;
            bins = B;

            N = Problem.N;
            D = Problem.D;

            Kmin = max(1, round(Kmin));
            Kmax = max(Kmin, round(D * min(1, max(0.01, KmaxFrac))));

            S = randn(N, D);
            K = randi([Kmin, Kmax], N, 1);
            Algorithm.freq = 0.5 * ones(1, D);

            Dec = Algorithm.DecodeTopK(S, K);
            Pop = Problem.Evaluation(Dec);

            while Algorithm.NotTerminated(Pop)
                Algorithm.freq = Algorithm.UpdateFreq(Pop, Algorithm.freq, rho);

                eliteIdx = Algorithm.SelectElites(Pop, memN);
                [MemDec, MemS, MemK] = Algorithm.MemeticSwap( ...
                    S(eliteIdx,:), K(eliteIdx), Algorithm.freq, Kmin, Kmax);

                if ~isempty(MemDec)
                    MemPop = Problem.Evaluation(MemDec);
                else
                    MemPop = [];
                end

                OffS = zeros(N, D);
                OffK = zeros(N, 1);
                bias = gamma * (2 * Algorithm.freq - 1);

                for i = 1:N
                    pool = [1:i-1, i+1:N];
                    r = pool(randperm(N-1, 3));
                    r1 = r(1);
                    r2 = r(2);
                    r3 = r(3);

                    vS = S(r1,:) + F * (S(r2,:) - S(r3,:)) + bias + sigma * randn(1, D);

                    vK = round(K(r1) + Fk * (K(r2) - K(r3)) + 0.5 * randn());
                    vK = min(Kmax, max(Kmin, vK));

                    uS = S(i,:);
                    jrand = randi(D);
                    mask = rand(1, D) < CR;
                    mask(jrand) = true;
                    uS(mask) = vS(mask);

                    if rand < 0.5
                        uK = vK;
                    else
                        uK = K(i);
                    end

                    uS = min(10, max(-10, uS));
                    OffS(i,:) = uS;
                    OffK(i)   = uK;
                end

                OffDec = Algorithm.DecodeTopK(OffS, OffK);
                OffPop = Problem.Evaluation(OffDec);

                if isempty(MemPop)
                    CandPop = [Pop, OffPop];
                    CandS   = [S; OffS];
                    CandK   = [K; OffK];
                else
                    CandPop = [Pop, OffPop, MemPop];
                    CandS   = [S; OffS; MemS];
                    CandK   = [K; OffK; MemK];
                end

                [Pop, S, K] = Algorithm.EnvSelectStratified( ...
                    CandPop, CandS, CandK, N, bins, Kmin, Kmax);
            end
        end
    end

    methods(Access = private)
        function Dec = DecodeTopK(~, S, K)
            [N, D] = size(S);
            Dec = false(N, D);

            for i = 1:N
                k = K(i);
                [~, idx] = sort(S(i,:), 'descend');
                Dec(i, idx(1:k)) = true;
            end
        end

        function fNew = UpdateFreq(~, Pop, fOld, rho)
            objs = Pop.objs;
            Cons = Pop.cons;
            hasCon = ~isempty(Cons);

            if hasCon
                CV = sum(max(0, Cons), 2);
                feas = find(CV == 0);
                if ~isempty(feas)
                    use = feas;
                else
                    [~, ord] = sort(CV, 'ascend');
                    use = ord(1:min(ceil(numel(ord) * 0.3), numel(ord)));
                end
                [FNo, ~] = NDSort(objs(use,:), Cons(use,:), numel(use));
            else
                use = 1:size(objs, 1);
                [FNo, ~] = NDSort(objs(use,:), numel(use));
            end

            eliteUse = use(FNo <= 2);
            if isempty(eliteUse)
                fNew = fOld;
                return;
            end

            Decs = Pop.decs;
            Decs = Decs(eliteUse,:);
            freqNow = mean(Decs > 0.5, 1);
            fNew = (1 - rho) * fOld + rho * freqNow;
            fNew = min(0.99, max(0.01, fNew));
        end

        function eliteIdx = SelectElites(~, Pop, m)
            objs = Pop.objs;
            Cons = Pop.cons;
            hasCon = ~isempty(Cons);

            if hasCon
                CV = sum(max(0, Cons), 2);
                feas = find(CV == 0);
                if ~isempty(feas)
                    [~, ord] = sort(objs(feas, 2), 'ascend');
                    eliteIdx = feas(ord(1:min(m, numel(feas))));
                    return;
                end
                [~, ord] = sort(CV, 'ascend');
                eliteIdx = ord(1:min(m, numel(ord)));
            else
                [~, ord] = sort(objs(:, 2), 'ascend');
                eliteIdx = ord(1:min(m, numel(ord)));
            end
        end

        function [MemDec, MemS, MemK] = MemeticSwap(~, S_elite, K_elite, freq, Kmin, Kmax)
            [m, D] = size(S_elite);
            MemDec = [];
            MemS   = [];
            MemK   = [];

            for i = 1:m
                k = min(Kmax, max(Kmin, K_elite(i)));
                [~, idx] = sort(S_elite(i,:), 'descend');
                sel   = idx(1:k);
                unsel = idx(k+1:end);

                if isempty(unsel)
                    continue;
                end

                rem = sel(end);
                [~, bestU] = max(freq(unsel));
                add = unsel(bestU);

                if add == rem
                    continue;
                end

                s2 = S_elite(i,:);
                s2(rem) = s2(rem) - 1.0;
                s2(add) = s2(add) + 1.0;

                dec2 = false(1, D);
                dec2(sel) = true;
                dec2(rem) = false;
                dec2(add) = true;

                MemDec = [MemDec; dec2]; %#ok<AGROW>
                MemS   = [MemS; s2];     %#ok<AGROW>
                MemK   = [MemK; k];      %#ok<AGROW>
            end
        end

        function [Pop, S, K] = EnvSelectStratified(~, CandPop, CandS, CandK, N, bins, Kmin, Kmax)
            PopObj = CandPop.objs;
            Cons = CandPop.cons;

            if ~isempty(Cons)
                [FrontNo, ~] = NDSort(PopObj, Cons, N);
            else
                [FrontNo, ~] = NDSort(PopObj, N);
            end

            CrowdDis = CrowdingDistance(PopObj, FrontNo);
            [~, ord] = sortrows([FrontNo(:), -CrowdDis(:)]);
            ord = ord(:);

            bins  = max(2, round(bins));
            edges = round(linspace(Kmin, Kmax + 1, bins + 1));
            cap   = ceil(N / bins);
            count = zeros(1, bins);

            chosen  = zeros(N, 1);
            pickedN = 0;

            for t = 1:numel(ord)
                idx = ord(t);
                k = CandK(idx);
                b = find(k >= edges(1:end-1) & k < edges(2:end), 1, 'first');

                if isempty(b)
                    b = min(bins, max(1, round((k - Kmin) / (Kmax - Kmin + eps) * bins) + 1));
                end

                if count(b) < cap
                    pickedN = pickedN + 1;
                    chosen(pickedN) = idx;
                    count(b) = count(b) + 1;
                    if pickedN >= N
                        break;
                    end
                end
            end

            if pickedN < N
                rest = setdiff(ord, chosen(1:pickedN), 'stable');
                chosen(pickedN+1:N) = rest(1:(N - pickedN));
            end

            Pop = CandPop(chosen);
            S   = CandS(chosen,:);
            K   = CandK(chosen);
        end
    end
end
