clear all; close all

activity_all=readtable('Circadian_time_series.xlsx');

name={};
treatment={};
v2=[];
for o=1:3
    
    ap={'light', 'LD', 'dark'};
    activity_all_regime = activity_all(:, contains(activity_all.Properties.VariableNames, char(ap(o))));

    name={};
    treatment={};
    circadianstrength=[];
    STACstrength=[];
    period9=[];
    for m=1:width(activity_all_regime)

        activity=activity_all_regime(:,m);
        name=vertcat(cellstr(name),cellstr(activity.Properties.VariableNames{1}));
        activity=table2array(activity); 

        nans(m)=sum(isnan(activity));
        infs(m)=sum(isinf(activity));

        f12=smoothdata(activity,'gaussian',15);
        f12=rescale(f12);
        
        power=[];
        [dpoaeCWT,f,coi] = cwt(f12, 1/30); % wavelet analysis
        for j=1:length(f12)
            cfsOAE = dpoaeCWT(:,j);
            q=abs(cfsOAE);
            h=horzcat(q,f);
            ex=coi(j);
            q2 = h(:,2) < ex; % exclude data points in the cone of influence 
            q(q2) = NaN;
            q3(:,j) = q;
        end

        for i=1:height(q3)
            power(i)=sum(q3(i,:), [], "omitnan");
        end
        wperiod=hours(seconds(1./f));
        [minDistance, indexOfMin] = min(abs(wperiod-24));
        circadianstrength(m)=power(indexOfMin);

        [M,Y] = max(power)
        period9(m)=1/f(Y)*0.01667
        STACstrength(m)=M;

        treatment=vertcat(treatment,char(ap(o)));

        %%%% Plotting %%%%

        % figure
        % subplot(1,2,1)
        % plot(f12)
        % subplot(1,2,2)
        % plot(hours(seconds(1./f)),power)

    end

    v=table(name,treatment, circadianstrength', period9', STACstrength');

    v2=vertcat(v2,v);

end

v2.Properties.VariableNames ={'species' 'treatment' 'circadian_strength' 'waveperiod' 'ultradian_strength'}






clear all; close all

activity_all=readtable('Circadian_time_series.xlsx');
activity_all_regime = activity_all(:, contains(activity_all.Properties.VariableNames, 'LD'));
    
name={};
treatment={};
v3=[];
for m=1:width(activity_all_regime)
 
activityt4=activity_all_regime(:,m);
[name{1:8}] = deal(cellstr(activityt4.Properties.VariableNames{1}));
activityt4=table2array(activityt4); 

time=(1:length(activityt4))*30*0.0002778;

f12=smoothdata(activityt4,'gaussian',15);
f12=rescale(f12);
 
phototime=horzcat(mean(time(1:1101)), mean(time(1102:2541)), mean(time(2542:3981)), mean(time(3982:5421)), mean(time(5422:6861)), mean(time(6862:8300)), mean(time(8301:9740)), mean(time(9741:11180)));
photoactive=horzcat(mean(f12(1:1101)), mean(f12(1102:2541)), mean(f12(2542:3981)), mean(f12(3982:5421)), mean(f12(5422:6861)), mean(f12(6862:8300)), mean(f12(8301:9740)), mean(f12(9741:11180)));

z9=[1,1101,1102,2541,2542,3981,3982,5421,5422,6861,6862,8300,8301,9740,9741,11180];
ppa=1;
ppb=2;
power=[];
for pp=1:length(phototime)
    fpp=f12(z9(ppa):z9(ppb));
    fpp=rescale(fpp);    
    [dpoaeCWT,f,coi] = cwt(fpp, 1/30);
    t=linspace(0,(length(fpp)*30*2.7778e-04),length(dpoaeCWT))';
    [coil, index] = unique(coi); 
    YIsx1 = interp1(coil, t(index), f);
    [coir, index] = unique(coi(((length(fpp)/2)):end)); 
    YIsx2 = interp1(coir, t(index)+(t(end)/2), f);
    power=[];
    q3=[];
        for j=1:length(fpp)
            cfsOAE = dpoaeCWT(:,j);
            q=abs(cfsOAE);
            h=horzcat(q,f);
            ex=coi(j);
            q2 = h(:,2) < ex; % exclude data points in the cone of influence 
            q(q2) = NaN;
            q3(:,j) = q;
        end

        for i=1:height(q3)
            power(i)=sum(q3(i,:), [], "omitnan");
        end
        [M,Y] = max(power)
        periodpp(pp)=1/f(Y)*0.01667
        waveletpp(pp)=M;

    ppa=ppa+2;
    ppb=ppb+2;
end


%%%% Plotting %%%%

% diel=1:length(activityt4); 
% diel(1:1101) = 1;
% diel(1102:2541) = 0;
% diel(2542:3981) = 1;
% diel(3982:5421) = 0;
% diel(5422:6861) = 1;
% diel(6862:8300) = 0;
% diel(8301:9740) = 1;
% diel(9741:11180) = 0;
% diel(11181:end) = 1;
% 
% f200=smoothdata(activityt4,'gaussian',1000);
% f200=rescale(f200);
% 
% figure
% plot(time,f12,'k','LineWidth',2)
% [pks1,locs1,w1,p1] = findpeaks(f12,'MinPeakProminence',.2);
% hold on
% area(time,diel, 'FaceColor', 'yellow', 'Facealpha', .3)
% hold on
% plot(time,f200,'g','LineWidth',3)
% hold on
% plot(phototime,photoactive,'-s' ,'LineWidth', 3, 'MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
% yyaxis right
% plot(phototime,waveletpp,'-s' ,'LineWidth', 3, 'MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor',[1 .6 .6])
% yyaxis left
% set(gca,'fontsize',15)
% set(gca,'linewidth',2)
% set(gca,'TickDir','out')
% set(gca, 'FontName', 'Times')
% a = gca;
% a.Color = 'k';
% set(a,'box','off','color','none')
% b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
% set(b,'linewidth',2)
% axes(a)
% linkaxes([a b])
% ylabel('Proportion of colony changed')
% title(char(string(name(1))))



TOD= {'Day' 'Night' 'Day' 'Night' 'Day' 'Night' 'Day' 'Night'};
v=table(string(name'),TOD', periodpp', waveletpp', photoactive');

v3=vertcat(v3,v);

end

v3.Properties.VariableNames ={'Colony' 'Time_of_day' 'periodpp' 'waveletpp' 'photoactive'}
