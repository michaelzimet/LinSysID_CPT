%PHYSIO System ID, after running the file PHYSIO_Plot_Main_Figs_Analyses.m

nx=3;
meas = 'BP'; %HR, CO, BP, SV

% step function input,
% prep time is 40s through 64s
%immersion interval 65s through 154s
% recovery is 155s to 190s
CPT_step=zeros(190,1);
CPT_step(65:154)=1;
ts_time = 65:190;
CPT_input = CPT_step(ts_time);

K = nx^2 + 2*nx;


loc = 'northeast';
switch meas
    case 'HR'
        mat = all_HR;
        thisYLabel='BPM (norm)';
        CPT_step=zeros(195,1);
        CPT_step(65:154)=1;
        ts_time = 65:195;
        CPT_input = CPT_step(ts_time);
    case 'CO'
        mat = all_CO;
        thisYLabel='L/min (norm)';
    case 'SV'
        mat = all_SV;
        thisYLabel='ml/m^2 (norm)';
    case 'BP'
        mat = all_BP;
        thisYLabel='mmHg (norm)';
        loc = 'southeast';
end

a = figure;
b = figure;
f = figure;
AICvec = zeros(5,2);
MSEvec = zeros(5,2);

for(iCond=1:1)
    switch iCond
        case 1
            condstr = 'Treatment';
        case 2
            condstr = 'Control';
    end

    steadvec=zeros(5,1);
    steadvec_disc=zeros(5,1);
    eigavec = zeros(nx,5);

    for(iOrder=1:1)
     % select line colors
            if      iOrder==1; thisColor = [255,0,0]; %red
            elseif  iOrder==2; thisColor = [255,140,0];% orange
            elseif  iOrder==3; thisColor = [252,226,5]; % yellow
            elseif  iOrder==4; thisColor = [0,255,0]; % green
            elseif  iOrder==5; thisColor = [0,0,255]; %blue
            end
            
  
        this_ts = squeeze(mean(mat(:, iOrder, iCond, ts_time), 'omitnan'));
        data = iddata(this_ts, CPT_input, 1);
        n4sid_opts = n4sidOptions('N4Weight','CVA', 'N4Horizon',[20], 'EnforceStability', false, 'Focus', 'simulation');
        sys = n4sid(data, nx,n4sid_opts,'DisturbanceModel','none');
        [ysys,fit,ic]=compare(data,sys);
       
        %Resid = this_ts - ysys.OutputData;
        %MSE = (1/length(this_ts))*sum(Resid' * Resid);
        MSE = sys.Report.Fit.MSE;
        AIC = length(this_ts)*log(MSE) + 2*K + length(this_ts)*(log(2*pi)+1);
        MSEvec(iOrder, iCond) = MSE;
        AICvec(iOrder, iCond) = AIC;
        

        %start at t=0
        figure(a)
        subplot(2, 1, iCond)
        plot(0:length(ysys.OutputData)-1, ysys.OutputData, 'Color', thisColor./255, 'LineWidth', 3)
        xlabel('Time (s)','Interpreter', 'latex')
        ylabel([meas ', ' thisYLabel], 'Interpreter','latex');
        ylim([0 1])
        xlim([0 130])
        xticks(0:30:120)
        yticks(0:0.2:1.0)
        hold on

        figure(b)
        subplot(5,2,iOrder*2 - (iCond ==1))
        plot(0:length(ysys.OutputData)-1,this_ts);
        hold on
        plot(0:length(ysys.OutputData)-1, ysys.OutputData);
        %xlabel('Time (s)','Interpreter', 'latex')
        %ylabel([meas], 'Interpreter','latex');
        ylim([0 1])
        xlim([0 130])
        xticks(0:30:120)
        yticks(0:0.5:1.0)
      
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
        
        %thisTransient = thisC * inv(thisA)*exp(20000*thisA)*thisB;
        %transientvec(iOrder) = thisTransient;

        %eigavec(:, iOrder)= eig(thisA); 

        %con_thisStead = con_sys.D - con_sys.C* inv(con_sys.A) * con_sys.B;
        %con_steadvec(iOrder)=con_thisStead;
        %con_eigavec(:, iOrder) = eig(con_sys.A);

        %Dynamical System Coefficients
        %[r1,r2,ru] = dyns(sys);
        %r1vec(iOrder) = r1;
        %r2vec(iOrder) = r2;
        %ruvec(iOrder) = ru;

        figure(f)
        subplot(2, 1, iCond)
        plot(0:89,ysys.OutputData(1:90), 'Color', thisColor./255, 'LineWidth',3)       
        xlabel('Time (s)','Interpreter', 'latex')
        %if(iCond==1)
        ylabel([meas ', ' thisYLabel], 'Interpreter','latex');
        %end
        ylim([0 1])
        yticks(0:0.2:1.0)
        xticks(0:30:90)
        hold on
        plot(100, thisStead, '.', "MarkerSize",40,'Color', thisColor./255)
        hold on
        end
      end
%AddLetters2Plots(a, {'A', 'B'},  'HShift', -0.08, 'VShift', -0.015)
saveas(a, ['figs/' meas '_models.eps'],'epsc')

%AddLetters2Plots(b, {'A', 'B','C','D','E','F','G','H', 'I', 'J'},  'HShift', -0.08, 'VShift', -0.015)
saveas(b, ['figs/' meas '_physioSysID.eps'],'epsc')

%AddLetters2Plots(f, {'A', 'B'},  'HShift', -0.08, 'VShift', -0.015)
saveas(f, ['figs/stead.eps'],'epsc')


