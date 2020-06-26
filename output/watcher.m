 function watcher()
     i=0;
     path(path, '/home/dmalt/Code/matlab/gopsiicos/gopsiicos_source/source');
     ISD = load('InputData4Simulations.mat');
     GainSVDTh = 0.01;
     HM = ups.sim.PrepareTestForward(GainSVDTh);
     Ctx = ISD.Ctx;
     [~, change_date] = system('stat -c %y A_reduced.mat');
     f1 = figure;
     f2 = figure;
     while i < 1000000
        pause(1);
        [~, change_date_new] = system('stat -c %y A_reduced.mat');
        if ~strcmp(change_date_new, change_date)
            load A_reduced;
            clear con_idx;
            for i_pair = 1:length(A_reduced)
                [i,j] = linToSq(A_reduced(i_pair), length(HM.GridLoc));
                con_idx(i_pair, :) = [i, j];
            end
            disp(length(con_idx))
            close(f1)
            con = ups.Bundles(con_idx, HM, Ctx);
            f1 = figure;
            con.Plot(0.05);
            material dull;
            close(f2)
            f2 = figure;
            plot(abs(X_a_output)')
            change_date = change_date_new;
        end
     end
 end
