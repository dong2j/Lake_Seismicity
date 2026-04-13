clc; clear; close all;

D=load('./Raw_Data/Global.txt');
Time=D(:,1);
Mag=D(:,5);
Lat=D(:,2);
Lon=D(:,3);
Depth=D(:,4);

time_cur=[1987 1990 2000 2010 2020];
figure('Units','normalized', 'Position', [0.1 0.1 0.6 0.2]);

for i=1:length(time_cur)-1
    

    jkf=find(Time>=time_cur(i)&Time<time_cur(i+1)&Mag>=(5-1e-12)&Depth<=70);

    Mag_sel = Mag(jkf);
    Lat_sel = Lat(jkf);
    Lon_sel = Lon(jkf);
    Time_sel = Time(jkf);


    % Magnitude processing
    M0=Mag_sel;
    T0=Time_sel;
    diffm = diff(M0);
    nonzero_elements = diffm(diffm ~= 0);
    delta = min(abs(nonzero_elements)) / 2;
    dm = 0.1;

    % MLE estimate
    m = floor(min(M0)/dm)*dm : dm : max(M0);
    n0 = hist(M0, m);
    [~, ind] = max(n0);
    Mc_MLE = m(ind) + 0.2;
    [Mc_MBS, b_es, number] = MBS_MLE_discrete(M0, delta);
    mc = max([Mc_MLE, Mc_MBS]);

    idx = M0 - mc >= -1e-11;

    [Dt, idx] = sort(T0(idx));   % Sort Dt and get sorting indices
    Dm = M0(idx) - mc;           % Apply the same order to Dm


    mmin = 0;
    [b_GR, loglikelihood_GR] = Estimation_GR_discrete(Dm, mmin, delta);

    subplot(1, length(time_cur)-1,i)

      plt_MFD(M0,Dm,mc,b_GR)
      fortitle=[num2str(time_cur(i)), ' to ',num2str(time_cur(i+1))];
      title(fortitle);
end
    
    
function []=plt_MFD(M,Dm,mc,b_GR)
    %% Plot cum- and non-culm distribution
    L=length(M);
    dm=0.1;
    m=floor(min(M)/dm)*dm:dm:max(M);
    if max(m)<max(M)
        m=[m,max(M)];
    end
    n0=hist(M,m);   
    for i=1:length(m)
        cn0(i)=length(find((M-m(i))>=-1e-10));
    end
    n=log10(n0);
    cn=log10(cn0);
    color1=[0.75 0.75 0.75; 0.25 0.25 0.25];
    color2=1/255*[132 94 194;178 91 0;0 139 200];
    plot_n0=semilogy(m,n0,'o', 'color', 0.4*[1 1 1], ...
                'markersize', 6.5, 'markerfacecolor', color1(1,:), ...
                'MarkerEdgeColor', color1(1,:));
    hold on;
    plot_cn0=semilogy(m,cn0,'o', 'color', 0.4*[1 1 1], ...
                'markersize', 6.5, 'markerfacecolor', color1(2,:), ...
                'MarkerEdgeColor', color1(2,:));
    hold on;
    
    if isnan(mc)==0
        %% Plot the best fit
        m=0:dm:max(Dm);
        if max(m)<max(Dm)
            m=[m,max(Dm)];
        end
        Moment_min=10^(1.5*(min(m)+6.07));
        F_GR=(Moment_min./10.^(3/2*(m+6.07))).^(2/3*b_GR);

        
        fit_GR=semilogy(m+mc,length(Dm)*F_GR,'color',color2(3,:),'linewidth',3);
        hold on;
        scatter(mc,length(Dm)*F_GR(1)*1.5,100,1/255*[255,165,0],'v','filled');
        hold on;
    
 
        %%
        % Set x-axis limits and ticks
        xMin = 5;
        xMax = 10;
        xlim([xMin xMax]);
        xticks(xMin:xMax); % Sets ticks at every integer between xMin and xMax
        xticklabels(arrayfun(@num2str, xMin:xMax, 'UniformOutput', false));

        
        % Set y-axis limits and ticks
        yMin = 0; % Assuming the minimum y value starts from 1
        yMax = log10(2e4);  % L is the maximum value in the original y data
        ylim([10^yMin, 10^yMax]);

        % Set y-axis ticks and labels
        yticks(10.^(yMin:yMax)); % Logarithmic ticks
        yticklabels(arrayfun(@(y) sprintf('10^{%d}', round(log10(y))), yticks, 'UniformOutput', false));

        
        ax = gca; 
        if ax.XLim(2)>0
            xTextPos = 0.95*ax.XLim(2); 
        else
            xTextPos = 1.05*ax.XLim(2); 
        end
        yTextPos = 0.9*10^(log10(ax.YLim(2)));

textStr = sprintf('$m_{\\mathrm{C}}$: %.1f, $b$: %.1f', mc, b_GR);
      
        text(xTextPos, yTextPos, textStr,'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 16,'Interpreter', 'latex');


    end
    grid on;box off;
        xlabel('\it m');
        ylabel('\it N');
    set(gca,'fontsize',13);
    hold on;   
end

function [Mc,b_es,number]=MBS_MLE_discrete(Dm,delta)
    dm=0.1;
    D00=sort(Dm,'descend');
    m=ceil(min(D00)/dm)*dm:dm:floor(D00(50)/dm)*dm;
    for k=1:length(m)
        Mmin=m(k);
        jkf=Dm>=(Mmin-1e-11);
        Dm0=Dm(jkf);
        b0(k)=bmle_discrete(Dm0,Mmin,delta);
        sigma(k)=b0(k)/sqrt(length(Dm0));
    end
    Delta=0.5;
    number=Delta/dm+1;
    b0_mean=[];
    for k=1:length(m)-number+1
        b0_mean(k)=mean(b0(k:k+number-1));
    end
    b02=b0(1:length(m)-number+1);
    threshold=abs(b02-b0_mean)-sigma(1:length(m)-number+1);
    jkf=find(threshold<0);
    if ~isempty(jkf)
        Mc=m(jkf(1));
        jkff=find(Dm>=(Mc-1e-11));
        number=length(jkff);
        b_es=b02(jkf(1));
    else
        [~,ind]=min(threshold);
        Mc=m(ind);
        jkff=find(Dm>=(Mc-1e-11));
        number=length(jkff);
        b_es=b02(ind);
    end
end

function [b_estimation_MLE]=bmle_discrete(Dm0,Mmin,delta)
    b_estimation_MLE=1/(log(10)*(mean(Dm0)-Mmin+delta));
end


function [b_GR,loglikelihood_GR]=Estimation_GR_discrete(Dm,Mmin,delta)
    DMoment=10.^(1.5*(Dm+6.07));
    Moment_min0=10^(1.5*(Mmin+6.07));
    Moment_min=10^(1.5*((Mmin-delta)+6.07));
    jkf=find((DMoment-Moment_min0)>=-1e-16);
    n=length(jkf);
    beta_GR=n/sum(log( DMoment(jkf)/Moment_min ));
    b_GR=beta_GR*3/2;
    loglikelihood_GR=n*(beta_GR*log(Moment_min)+log(beta_GR))-(1+beta_GR)*sum(log(DMoment(jkf)));
end
