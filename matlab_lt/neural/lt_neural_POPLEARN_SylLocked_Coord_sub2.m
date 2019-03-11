            %%
            figcount=1;
            subplotrows=5;
            subplotcols=2;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            [~, indsort] = sort(cohscal_tmp);
            nplot = 4;
            
            % --- collect all rho and power to plot in end
%             rhotmp = [];
%             powtmp = [];
%             frtmp = [];
            cohtmp = [];
            
            % ====== LOW CORREALTION TRIALS
            trialstoplot = indsort(1:nplot);
            ptit = 'LOW COHERNECE';

            for tt=trialstoplot
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots; hsplot];
                title([ptit ', tr#' num2str(tt)]);
                ylabel('0.1*lfp; spikes');
                
                coh = cohscal_tmp(tt);
                
                % ================ LMAN
                pcol = [0.2 0.6 0.2];
                spk = cellfun(@(x)(x(tt)), spkdat_LMAN);
                lfp = cellfun(@(x)x(:, tt), lfpdat_LMAN,  'UniformOutput', 0);
                
                % -- plot spikes
                for j=1:length(spk)
                    lt_neural_PLOT_rasterline(spk{j}, j, pcol);
                end
                % -- plot lfp
                xlfp = PARAMS.THIS.lfpx;
                for j=1:length(lfp)
                    y = 0.1*lfp{j};
                    plot(xlfp, y, 'Color', pcol);
                end
                
                % ================ RA
                pcol = 'r';
                spk = cellfun(@(x)(x(tt)), spkdat_RA);
                lfp = cellfun(@(x)x(:, tt), lfpdat_RA,  'UniformOutput', 0);
                yshift = 5;
                
                % -- plot spikes
                for j=1:length(spk)
                    lt_neural_PLOT_rasterline(spk{j}, j+yshift, pcol);
                end
                % -- plot lfp
                xlfp = PARAMS.THIS.lfpx;
                for j=1:length(lfp)
                    y = 0.1*lfp{j};
                    plot(xlfp, y+yshift, 'Color', pcol);
                end
                
                % -- annotate rho and power
                tstr = ['coh = ' num2str(coh)];
                lt_plot_annotation(1, tstr, 'b');
                cohtmp = [cohtmp; coh];
                lt_plot_zeroline_vert;
                lt_plot_zeroline;
            end
            
            % ====== MEDIAN CORREALTION TRIALS
            trialstoplot = indsort(floor(length(indsort)/2):floor(length(indsort)/2)+nplot-3);
            ptit = 'MID COHERNECE';
            
            for tt=trialstoplot
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots; hsplot];
                title([ptit ', tr#' num2str(tt)]);
                ylabel('0.1*lfp; spikes');
                
                coh = cohscal_tmp(tt);
                
                % ================ LMAN
                pcol = [0.2 0.6 0.2];
                spk = cellfun(@(x)(x(tt)), spkdat_LMAN);
                lfp = cellfun(@(x)x(:, tt), lfpdat_LMAN,  'UniformOutput', 0);
                
                % -- plot spikes
                for j=1:length(spk)
                    lt_neural_PLOT_rasterline(spk{j}, j, pcol);
                end
                % -- plot lfp
                xlfp = PARAMS.THIS.lfpx;
                for j=1:length(lfp)
                    y = 0.1*lfp{j};
                    plot(xlfp, y, 'Color', pcol);
                end
                
                % ================ RA
                pcol = 'r';
                spk = cellfun(@(x)(x(tt)), spkdat_RA);
                lfp = cellfun(@(x)x(:, tt), lfpdat_RA,  'UniformOutput', 0);
                yshift = 5;
                
                % -- plot spikes
                for j=1:length(spk)
                    lt_neural_PLOT_rasterline(spk{j}, j+yshift, pcol);
                end
                % -- plot lfp
                xlfp = PARAMS.THIS.lfpx;
                for j=1:length(lfp)
                    y = 0.1*lfp{j};
                    plot(xlfp, y+yshift, 'Color', pcol);
                end
                
                % -- annotate rho and power
                tstr = ['coh = ' num2str(coh)];
                lt_plot_annotation(1, tstr, 'b');
                cohtmp = [cohtmp; coh];
                lt_plot_zeroline_vert;
                lt_plot_zeroline;
            end
            
            
            
            
                        % ====== HIGH CORREALTION TRIALS
            trialstoplot = indsort(end-nplot+1:end);
            ptit = 'HIGH COHERNECE';
            
            for tt=trialstoplot
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots; hsplot];
                title([ptit ', tr#' num2str(tt)]);
                ylabel('0.1*lfp; spikes');
                
                coh = cohscal_tmp(tt);
                
                % ================ LMAN
                pcol = [0.2 0.6 0.2];
                spk = cellfun(@(x)(x(tt)), spkdat_LMAN);
                lfp = cellfun(@(x)x(:, tt), lfpdat_LMAN,  'UniformOutput', 0);
                
                % -- plot spikes
                for j=1:length(spk)
                    lt_neural_PLOT_rasterline(spk{j}, j, pcol);
                end
                % -- plot lfp
                xlfp = PARAMS.THIS.lfpx;
                for j=1:length(lfp)
                    y = 0.1*lfp{j};
                    plot(xlfp, y, 'Color', pcol);
                end
                
                % ================ RA
                pcol = 'r';
                spk = cellfun(@(x)(x(tt)), spkdat_RA);
                lfp = cellfun(@(x)x(:, tt), lfpdat_RA,  'UniformOutput', 0);
                yshift = 5;
                
                % -- plot spikes
                for j=1:length(spk)
                    lt_neural_PLOT_rasterline(spk{j}, j+yshift, pcol);
                end
                % -- plot lfp
                xlfp = PARAMS.THIS.lfpx;
                for j=1:length(lfp)
                    y = 0.1*lfp{j};
                    plot(xlfp, y+yshift, 'Color', pcol);
                end
                
                % -- annotate rho and power
                tstr = ['coh = ' num2str(coh)];
                lt_plot_annotation(1, tstr, 'b');
                cohtmp = [cohtmp; coh];
                lt_plot_zeroline_vert;
                lt_plot_zeroline;
            end
            
            
            
            
            % ================ FORMAT ALL
            linkaxes(hsplots, 'xy');