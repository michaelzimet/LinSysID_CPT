%EYE System ID, after running the file EYE_Plot_Main_Figs_Analyses.m
nx=4;

% step function input,
% prep time is 40s through 64s
%immersion interval 65s through 154s
% recovery is 155s to 190s
CPT_step=zeros(195,1);
CPT_step(65:154)=1;

%modeling choices around prep and recovery
ts_time = 65:195;
%ts_time = 65:154;
CPT_input = CPT_step(ts_time);
K = nx^2 + 2*nx;

measures = {'Pupil'};
loc = 'northeast';
meas = measures{1};
mat = paMatAll;
thisYLabel='(norm)';

view = 'ss'; %'ss', 'ds', or 'tf'

a = figure;
t= tiledlayout(2,1);
aa(1) = nexttile;
title(aa(1),'Treatment','Interpreter','latex','FontSize',12)
hold(aa(1),'on')
aa(2) = nexttile;
title(aa(2),'Control','Interpreter','latex','FontSize',12)
hold(aa(2),'on')

b = figure;
tt= tiledlayout(5,2, TileIndexing="columnmajor");

%{
c = figure;
q= tiledlayout(2,1);
cc(1) = nexttile;
title(cc(1),'Treatment, 0.5 Hz','Interpreter','latex','FontSize',12)
hold(cc(1),'on')
cc(2) = nexttile;
title(cc(2),'Control, 0.5 Hz','Interpreter','latex','FontSize',12)
hold(cc(2),'on')
%}

d = figure;
qq= tiledlayout(5,2, TileIndexing="columnmajor");

e = figure;
r= tiledlayout(2,1);
ee(1) = nexttile;
title(ee(1),'Treatment','Interpreter','latex','FontSize',12)
hold(ee(1),'on')
ee(2) = nexttile;
title(ee(2),'Control','Interpreter','latex','FontSize',12)
hold(ee(2),'on')

f = figure;
ttt= tiledlayout(2,2);

AICvec = zeros(5,2);
MSEvec = zeros(5,2);
Resid_mean_vec = zeros(5,2);
Resid_std_vec = zeros(5,2);
Resid_mean_vec_ds = zeros(5,2);
Resid_std_vec_ds = zeros(5,2);


for(iCond=1:2)
    switch iCond
        case 1
            condstr = 'Treatment';
        case 2
            condstr = 'Control';
    end


    steadvec=zeros(5,1);
    transientvec=zeros(5,1);
    eigavec = zeros(nx,5);


   for(iOrder=1:5)
   % select line colors
        if      iOrder==1; thisColor = [255,0,0]; %red
        elseif  iOrder==2; thisColor = [255,140,0];% orange
        elseif  iOrder==3; thisColor = [252,226,5]; % yellow
        elseif  iOrder==4; thisColor = [0,255,0]; % green
        elseif  iOrder==5; thisColor = [0,0,255]; %blue
        end
                
        %average pupils over both eyes
        %average full population

        switch nx
            case 2
                horiz = [15, 20, 20];
                str = 'planar';
            otherwise
                horiz = [20];
                str = 'physio';
        end
        
        this_ts = squeeze(mean(mean(paMatAll(:, :, iOrder,iCond,ts_time),2), 'omitnan'));
        data = iddata(this_ts, CPT_input, 1);
        n4sid_opts = n4sidOptions('N4Weight','CVA', 'N4Horizon',horiz,'EnforceStability', false, 'Focus', 'simulation');
        sys = n4sid(data, nx,n4sid_opts,'DisturbanceModel','none');
        [ysys,fit,ic]=compare(iddata(this_ts, CPT_input,1),sys);

     
        Resid = this_ts - ysys.OutputData;
        MSE_calc = (1/length(this_ts))*sum(Resid' * Resid);
        MSE = sys.Report.Fit.MSE;
        AIC = length(this_ts)*log(MSE) + 2*K + length(this_ts)*(log(2*pi)+1);
        MSEvec(iOrder, iCond) = MSE;
        AICvec(iOrder, iCond) = AIC;

       %downsampling analysis
        this_ts_ds = this_ts(1:2:end);
        CPT_input_ds = CPT_input(1:2:end);
        data_ds = iddata(this_ts_ds, CPT_input_ds, 1);
        sys_ds = n4sid(data_ds, nx,n4sid_opts,'DisturbanceModel','none');
        [ysys_ds,fit_ds,ic_ds]=compare(data_ds,sys_ds);
       
        Resid_ds = this_ts_ds - ysys_ds.OutputData;
        %MSE = (1/length(this_ts))*sum(Resid' * Resid);
        %MSE_ds = sys_ds.Report.Fit.MSE;
        %AIC_ds = length(this_ts_ds)*log(MSE_ds) + 2*K + length(this_ts_ds)*(log(2*pi)+1);
        %MSEvec_ds(iOrder, iCond) = MSE_ds;
        %AICvec_ds(iOrder, iCond) = AIC_ds;
        %}
        
         %Residual mean and standard deviation
        Resid_mean_vec(iOrder, iCond) = mean(Resid);
        Resid_std_vec(iOrder, iCond) = std(Resid);
        Resid_mean_vec_ds(iOrder, iCond) = mean(Resid_ds);
        Resid_std_vec_ds(iOrder, iCond) = std(Resid_ds);


        %start at t=0
        figure(a)
        plot(aa(iCond),0:length(ysys.OutputData)-1, ysys.OutputData, 'Color', thisColor./255, 'LineWidth', 3)
        ylim(aa(1),[0 1])
        xlim(aa(1),[0 130])
        xticks(aa(1),0:30:120)
        yticks(aa(1),0:0.2:1.0)
        ylim(aa(2),[0 1])
        xlim(aa(2),[0 130])
        xticks(aa(2),0:30:120)
        yticks(aa(2),0:0.2:1.0)
        hold on
      
        %{
        figure(c)
        plot(cc(iCond),2*(0:length(ysys_ds.OutputData)-1), ysys_ds.OutputData, 'Color', thisColor./255, 'LineWidth', 3)
        ylim(cc(1),[0 1])
        xlim(cc(1),[0 130])
        xticks(cc(1),0:30:120)
        yticks(cc(1),0:0.2:1.0)
        ylim(cc(2),[0 1])
        xlim(cc(2),[0 130])
        xticks(cc(2),0:30:120)
        yticks(cc(2),0:0.2:1.0)
        hold on
        %}

        thisStead_disc = -sys.C * inv(sys.A - eye(nx)) * sys.B;
        steadvec_disc(iOrder) = thisStead_disc;

        %Continuous Time
        cont_sys = d2c(sys);
        [ysys_cont,fit_cont,ic_cont]=compare(data,cont_sys);

        thisA = cont_sys.A;
        thisB = cont_sys.B;
        thisC = cont_sys.C;
        thisD = cont_sys.D;

        thisStead = thisD - thisC* inv(thisA) * thisB;
        steadvec(iOrder) = thisStead;

        %eigavec(:, iOrder)= eig(thisA); 

        figure(b)
        ab(iOrder*2 - sum(iCond==1)) = nexttile;
        plot(0:length(ysys.OutputData)-1,this_ts)
        hold on
        plot(0:length(ysys.OutputData)-1, ysys.OutputData)
        ylim([0 1])
        xlim([0 130])
        xticks(0:30:120)
        yticks(0:0.5:1.0)

        
        figure(d)
        ad(iOrder*2 - sum(iCond==1)) = nexttile;
        plot(2*(0:length(ysys_ds.OutputData)-1),this_ts_ds)
        hold on
        plot(2*(0:length(ysys_ds.OutputData)-1), ysys_ds.OutputData)
        ylim([0 1])
        xlim([0 130])
        xticks(0:30:120)
        yticks(0:0.5:1.0)
        

        figure(e)
        %plot(ee(iCond),0:length(ysys.OutputData)-1, Resid, 'Color', thisColor./255, 'LineWidth', 3)
        %plot(ee(iCond),Resid_auto, 'Color', thisColor./255, 'LineWidth', 3)
        histogram(ee(iCond), Resid, 'BinWidth', 0.01, 'EdgeColor', thisColor./255,'FaceColor', thisColor./255,'FaceAlpha',0.05,  'LineWidth', 3, 'Normalization', 'percentage')
        ylim(ee(1),[0 35])
        xlim(ee(1),[-0.3 0.3])
        xticks(ee(1),-0.3:0.1:0.3)
        %yt = yticks;
        %yticklabels(ee(1), num2str(yt*100))
        %yticks(ee(1),0:0.2:1.0)
        ylim(ee(2),[0 35])
        xlim(ee(2),[-0.3 0.3])
        xticks(ee(2),-0.3:0.1:0.3)
        %yticklabels(ee(2), num2str(yticks*100))
        %yticks(ee(2),0:0.2:1.0)
        hold on



        figure(f)
        if iOrder ==1 
            nexttile;
        end
        if iCond == 1
            ylabel([meas ', ' thisYLabel], 'Interpreter','latex');
        end
        plot(0:89,ysys.OutputData(1:90), 'Color', thisColor./255, 'LineWidth',3)       
        ylim([0 1])
        yticks(0:0.2:1.0)
        xticks(0:30:90)
        hold on
        plot(100, thisStead, '.', "MarkerSize",40,'Color', thisColor./255)
        hold on
        end
      end


figure(a)
xlabel(t, 'Time (s)','Interpreter', 'latex')
ylabel(t,[meas ', ' thisYLabel], 'Interpreter','latex');
saveas(a, ['figs/' meas '_models.eps'],'epsc')

%{
figure(c)
xlabel(q, 'Time (s)','Interpreter', 'latex')
ylabel(q,[meas ', ' thisYLabel], 'Interpreter','latex');
saveas(c, ['figs/' meas '_models_ds.eps'],'epsc')
%}

figure(b)
for iOrder = 1:5
    ab_row = nexttile(iOrder);      % left tile of row
    ylabel(ab_row, "Trial " + num2str(iOrder), 'Interpreter','latex', 'FontSize',12)
end
title(ab(1),'Treatment','Interpreter','latex','FontSize',12)
title(ab(2),'Control','Interpreter','latex','FontSize',12)
xlabel(tt, 'Time (s)','Interpreter', 'latex')
ylabel(tt, [meas ', ' thisYLabel], 'Interpreter','latex')
saveas(b, ['figs/' meas '_physioSysID.eps'],'epsc')


figure(d)
for iOrder = 1:5
    ad_row = nexttile(iOrder);      % left tile of row
    ylabel(ad_row, "Trial " + num2str(iOrder), 'Interpreter','latex', 'FontSize',12)
end
title(ad(1),'Treatment, 0.5 Hz','Interpreter','latex','FontSize',12)
title(ad(2),'Control, 0.5 Hz','Interpreter','latex','FontSize',12)
xlabel(qq, 'Time (s)','Interpreter', 'latex')
ylabel(qq, [meas ', ' thisYLabel], 'Interpreter','latex')
saveas(d, ['figs/' meas '_physioSysID_ds.eps'],'epsc')


figure(e)
xlabel(r, 'Residual','Interpreter', 'latex')
ylabel(r, 'Frequency (\%)', 'Interpreter','latex');
saveas(e, ['figs/' meas '_Resid.eps'],'epsc')

figure(f)
xlabel(ttt, 'Time (s)','Interpreter', 'latex')
%saveas(f, ['figs/stead.eps'],'epsc')