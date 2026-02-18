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
b = figure;
tt= tiledlayout(5,2);
f = figure;
ttt= tiledlayout(2,2);
AICvec = zeros(5,2);
MSEvec = zeros(5,2);

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

     
        %Resid = this_ts - ysys.OutputData;
        %MSE = (1/length(this_ts))*sum(Resid' * Resid);
        MSE = sys.Report.Fit.MSE;
        AIC = length(this_ts)*log(MSE) + 2*K + length(this_ts)*(log(2*pi)+1);
        MSEvec(iOrder, iCond) = MSE;
        AICvec(iOrder, iCond) = AIC;


        %start at t=0
        figure(a)
        if iOrder==1
            nexttile;
        end
        plot(0:length(ysys.OutputData)-1, ysys.OutputData, 'Color', thisColor./255, 'LineWidth', 3)
        ylim([0 1])
        xlim([0 130])
        xticks(0:30:120)
        yticks(0:0.2:1.0)
        hold on
      
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
        nexttile;      
        plot(0:length(ysys.OutputData)-1,this_ts)
        hold on
        plot(0:length(ysys.OutputData)-1, ysys.OutputData)
        ylim([0 1])
        xlim([0 130])
        xticks(0:30:120)
        yticks(0:0.5:1.0)

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

figure(b)
xlabel(tt, 'Time (s)','Interpreter', 'latex')
ylabel(tt, [meas ', ' thisYLabel], 'Interpreter','latex')
saveas(b, ['figs/' meas '_physioSysID.eps'],'epsc')

figure(f)
xlabel(ttt, 'Time (s)','Interpreter', 'latex')
%saveas(f, ['figs/stead.eps'],'epsc')
