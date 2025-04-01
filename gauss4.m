function [areacoverage, fa] = gauss4(time,signal,n,weights,cor_chrom_1, cor_chrom_2,cor_chrom_3, deconmode)
% [y_fit,yhat,areacoverage,fa]
% y_fit = total fit
% yhat = deconvoluted fits
% areacoverage = overlap fit vs measured distributin
% fa = fitted parameters (height, position, sigma)


[hgtt, loct, wdtt] = findpeaks(signal,time);%,'MinPeakHeight', (max(signal) - min(signal)) * 0.01 ,'NPeaks', n);
[~,perm] = sort(hgtt);
strT = 12.5;
ndT = 14.3;
hgtt = hgtt(perm(end-n+1:end));
loct = loct(perm(end-n+1:end));
wdtt = wdtt(perm(end-n+1:end));
p_var_2 = [-0.00352544207707122	0.0976618734670780	-0.668567622117514];    % system variance polynomial

% With weight factor
if strcmp(deconmode, 'gauss')
    if weights == true
    fo = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[zeros(1,n)';ones(1,n)'.*strT;zeros(1,n)'],... % minimum time point
        'Upper',[ones(1,n)'.*max(signal);ones(1,n)'.*ndT;ones(1,n)'.*0.1],...% max time point, max sigma (0.1 default)
        'StartPoint',[hgtt; loct; wdtt./2.35],...
        'Weights',signal.^2); % y^2 used as weight factor to put more emphasis on peak top

    % No weight factor (in case of difficult fit)
    else
    fo = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[zeros(1,n)';ones(1,n)'.*strT;zeros(1,n)'],... % minimum time point
        'Upper',[ones(1,n)'.*max(signal);ones(1,n)'.*ndT;ones(1,n)'.*0.1],...% max time point, max sigma
        'StartPoint',[hgtt; loct; wdtt./2.35]); 
    end

    x = time;
    if n == 1
        ft1 = fittype('(a1 / (c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b1) .^2) / (c1 .^ 2))','options', fo); % includes weight factor (in case if true)
        ft = ft1;
    elseif n ==2
        ft2 = fittype('(a1 / (c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b1) .^2) / (c1 .^ 2)) + (a2 / (c2 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b2) .^2) / (c2 .^ 2))','options', fo);
        ft = ft2;
    elseif n ==3
        ft3 = fittype('(a1 / (c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b1) .^2) / (c1 .^ 2)) + (a2 / (c2 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b2) .^2) / (c2 .^ 2)) + (a3 / (c3 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b3) .^2) / (c3 .^ 2))','options', fo);
        ft = ft3;
    elseif n == 4
        ft4 = fittype('(a1 / (c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b1) .^2) / (c1 .^ 2)) + (a2 / (c2 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b2) .^2) / (c2 .^ 2)) + (a3 / (c3 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b3) .^2) / (c3 .^ 2)) + (a4 / (c4 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b4) .^2) / (c4 .^ 2))','options', fo);
        ft = ft4;
    elseif n == 5
        ft5 = fittype('(a1 / (c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b1) .^2) / (c1 .^ 2)) +(a2 / (c2 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b2) .^2) / (c2 .^ 2)) + (a3 / (c3 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b3) .^2) / (c3.^ 2)) +(a4 / (c4 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b4) .^2) / (c4 .^ 2))+(a5 / (c5 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b5) .^2) / (c5 .^ 2))','options', fo);
        ft = ft5;
    else
        error("ERROR: too many peaks for fitting, please increase treshold")
    end

    fa =fit(time,signal,ft);

elseif or(strcmp(deconmode, 'stacked'), strcmp(deconmode, 'gausseqv'))
 
   fo = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[zeros(1,n)';ones(1,n)'.*strT;zeros(1,n)';ones(n,1).*1;ones(n,1).*-0.15],... % minimum time point
        'Upper',[ones(1,n)'.*max(signal);ones(1,n)'.*ndT;ones(1,n)'.*0.5;ones(n,1).*1000;ones(n,1).*0.15],...% max time point, max sigma
        'StartPoint',[hgtt; loct; wdtt./2.35; ones(n,1).*1000; zeros(n,1)]); 
    
    x = time;
    if n == 1
        ft1 = fittype('(a1 / (c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b1) .^2) / (c1 .^ 2))','options', fo); % includes weight factor (in case if true)
        sln1 = fittype('a1 .* (1 + ((x - b1) .^ 2 ./ (d1 .* (c1 + e1 .* (x - b1)) .^ 2))) .^ (-1 .* d1)', 'options', fo);
        ft = sln1;
        fg = ft1;
    elseif n ==2
        ft2 = fittype('(a1 / (c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b1) .^2) / (c1 .^ 2)) + (a2 / (c2 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b2) .^2) / (c2 .^ 2))','options', fo);
        sln2 = fittype('a1 .* (1 + ((x - b1) .^ 2 ./ (d1 .* (c1 + e1 .* (x - b1)) .^ 2))) .^ (-1 .* d1) + a2 .* (1 + ((x - b2) .^ 2 ./ (d2 .* (c2 + e2 .* (x - b2)) .^ 2))) .^ (-1 .* d2)', 'options', fo);
        ft = sln2;
        fg = ft2;
    elseif n ==3
        ft3 = fittype('(a1 / (c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b1) .^2) / (c1 .^ 2)) + (a2 / (c2 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b2) .^2) / (c2 .^ 2)) + (a3 / (c3 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b3) .^2) / (c3 .^ 2))','options', fo);
        sln3 = fittype('a1 .* (1 + ((x - b1) .^ 2 ./ (d1 .* (c1 + e1 .* (x - b1)) .^ 2))) .^ (-1 .* d1) + a2 .* (1 + ((x - b2) .^ 2 ./ (d2 .* (c2 + e2 .* (x - b2)) .^ 2))) .^ (-1 .* d2) + a3 .* (1 + ((x - b3) .^ 2 ./ (d3 .* (c3 + e3 .* (x - b3)) .^ 2))) .^ (-1 .* d3)', 'options', fo);
        ft = sln3;
        fg = ft3;
    elseif n == 4
        ft4 = fittype('(a1 / (c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b1) .^2) / (c1 .^ 2)) + (a2 / (c2 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b2) .^2) / (c2 .^ 2)) + (a3 / (c3 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b3) .^2) / (c3 .^ 2)) + (a4 / (c4 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b4) .^2) / (c4 .^ 2))','options', fo);
        sln4 = fittype('a1 .* (1 + ((x - b1) .^ 2 ./ (d1 .* (c1 + e1 .* (x - b1)) .^ 2))) .^ (-1 .* d1) + a2 .* (1 + ((x - b2) .^ 2 ./ (d2 .* (c2 + e2 .* (x - b2)) .^ 2))) .^ (-1 .* d2) + a3 .* (1 + ((x - b3) .^ 2 ./ (d3 .* (c3 + e3 .* (x - b3)) .^ 2))) .^ (-1 .* d3) + a4 .* (1 + ((x - b4) .^ 2 ./ (d4 .* (c4 + e4 .* (x - b4)) .^ 2))) .^ (-1 .* d4)', 'options', fo);
        ft = sln4;
        fg = ft4;
    else
        error("ERROR: too many peaks for fitting, please increase treshold")
    end

    fa =fit(time,signal,ft);
    slnpara = coeffvalues(fa);
    slnpara = reshape(slnpara,[],5);
    fas = [];
    for i = 1:n
        fog = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[zeros(1,1)';ones(1,1)'.*strT;zeros(1,1)'],... % minimum time point
            'Upper',[ones(1,1)'.*max(signal);ones(1,1)'.*ndT;ones(1,1)'.*0.1],...% max time point, max sigma
            'StartPoint',[slnpara(1,1);slnpara(1,2);slnpara(1,3)]); 
        gs1 = fittype('(a1 / (c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - b1) .^2) / (c1 .^ 2))','options', fog);
   
        y = SLN(time,slnpara(i,1),slnpara(i,2),slnpara(i,3),slnpara(i,4),slnpara(i,5));
        fag = fit(time,y, gs1);
        fas = [fas, coeffvalues(fag)];
    end
%     fax = reshape(fas,n,3)'
    
    if strcmp(deconmode, 'gausseqv')
        fa = [slnpara(:,1).*((slnpara(:,3)./sqrt(2)).* sqrt(2*pi())),slnpara(:,2),slnpara(:,3)./sqrt(2)];
    elseif strcmp(deconmode, 'stacked')
        fa = reshape(fas,n,3)';
    end
    
    fa = reshape(fa,1,[]);
   if n == 1 
       fa = cfit(fg,fa(1),fa(2),fa(3));
   elseif n == 2 
       fa = cfit(fg,fa(1),fa(2),fa(3),fa(4),fa(5),fa(6));
   elseif n == 3
       fa = cfit(fg,fa(1),fa(2),fa(3),fa(4),fa(5),fa(6),fa(7),fa(8),fa(9));    
   elseif n == 4
       fa = cfit(fg,fa(1),fa(2),fa(3),fa(4),fa(5),fa(6),fa(7),fa(8),fa(9),fa(10),fa(11),fa(12));    
   else
       error("ERROR: too many peaks for fitting, please increase treshold")
    end 
end
    
if n == 1
    ft = ft1;
    yhat = (fa.a1 / (fa.c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b1) .^2) / (fa.c1 .^ 2));
    y_fit = sum(yhat,2);
    figure
    hold on
    plot(x,signal, 'LineWidth', 2,'color','#0D2747');
    %plot(x,y_fit, 'LineWidth', 2,'color','#D06D91');       % fitted data
    plot(x,yhat(:,1), 'LineWidth', 2,'color','#567964','LineStyle','--');    % population 1
    plot(time,cor_chrom_1,'LineWidth',2,'color', '#567964'); % corrected chrom
    xlabel('Time (min)'); ylabel('Normalized intensity'); 
    set(gca, 'FontSize', 14, 'TickDir','out','LineWidth',2);
    ax = gca; 
    xlim([12.5,14.3]);
    ylim([0,1]);
    %legend('PSNP blend', 'Fitted curve','Population 1','PSD population 1');
    %legend ('boxoff');
    %legend ('location','northwest');
    hold off
    [fit_h, fit_loc, fit_FWHM] = findpeaks(y_fit,x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend');
    fit_para = [fit_h,fit_loc,fit_FWHM/2.35, polyval(p_var_2,fit_loc)]
    [p1_h, p1_loc, p1_FWHM] = findpeaks(yhat(:,1),x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend');
    p1_para = [p1_h, p1_loc, p1_FWHM/2.35, polyval(p_var_2,p1_loc)]
elseif n ==2
    ft = ft2;
    yhat1 = (fa.a1 / (fa.c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b1) .^2) / (fa.c1 .^ 2));
    yhat2 = (fa.a2 / (fa.c2 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b2) .^2) / (fa.c2 .^ 2));
    yhat =[yhat1 yhat2];
    y_fit = sum(yhat,2);
    figure
    hold on
    plot(x,signal, 'LineWidth', 2,'color','#0D2747');
    %plot(x,y_fit, 'LineWidth', 2,'color','#D06D91');       % fitted data
    plot(x,yhat(:,2), 'LineWidth', 2,'color','#567964','LineStyle','--');    % population 1
    plot(time,cor_chrom_1,'LineWidth',2,'color', '#567964'); % corrected chrom 1
    plot(x,yhat(:,1), 'LineWidth', 2,'color','#A7A7A7','LineStyle','--');    % population 2
    plot(time,cor_chrom_2,'LineWidth',2,'color', '#A7A7A7'); % corrected chrom 2
    hold off
    xlabel('Time (min)'); ylabel('Normalized intensity'); 
    set(gca, 'FontSize', 14, 'TickDir','out','LineWidth',2);
    ax = gca; 
    xlim([12.5,14.3]);
    ylim([0,1]);
    %legend('PSNP blend', 'Fitted curve','Population 1','PSD population 1','Population 2','PSD population 2');
    %legend('PSNP blend', 'Fitted curve','Population 1','PUR population 1','Population 2','PUR population 2');
    %legend ('boxoff');
    %legend ('location','northwest');
    [fit_h, fit_loc, fit_FWHM] = findpeaks(y_fit,x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend','NPeaks',2);
    fit_para = [fit_h,fit_loc,fit_FWHM/2.35, polyval(p_var_2,fit_loc)]
    [p1_h, p1_loc, p1_FWHM] = findpeaks(yhat(:,1),x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend');
    [p2_h, p2_loc, p2_FWHM] = findpeaks(yhat(:,2),x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend');
    p1_para = [p1_h, p1_loc, p1_FWHM/2.35, polyval(p_var_2,p1_loc)]
    p2_para = [p2_h, p2_loc, p2_FWHM/2.35, polyval(p_var_2,p2_loc)]

elseif n ==3
    ft = ft3;
    yhat1 = (fa.a1 / (fa.c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b1) .^2) / (fa.c1 .^ 2));
    yhat2 = (fa.a2 / (fa.c2 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b2) .^2) / (fa.c2 .^ 2));
    yhat3 = (fa.a3 / (fa.c3 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b3) .^2) / (fa.c3 .^ 2));
    yhat =[yhat1 yhat2 yhat3];
    y_fit = sum(yhat,2);
    figure
    hold on
    plot(x,signal, 'LineWidth', 2,'color','#0D2747');
    %plot(x,y_fit, 'LineWidth', 2,'color','#D06D91');       % fitted data
    plot(x,yhat(:,2), 'LineWidth', 2,'color','#567964','LineStyle','--');    % population 1
    plot(time,cor_chrom_1,'LineWidth',2,'color', '#A7A7A7'); % corrected chrom 1
    plot(x,yhat(:,1), 'LineWidth', 2,'color','#A7A7A7','LineStyle','--');    % population 2
    plot(time,cor_chrom_2,'LineWidth',2,'color', '#567964'); % corrected chrom 2
    plot(x,yhat(:,3), 'LineWidth', 2,'color','#9877AB','LineStyle','--');  % population 3
    plot(time,cor_chrom_3,'LineWidth',2,'color', '#9877AB'); % corrected chrom 3
    hold off
    xlabel('Time (min)'); ylabel('Normalized intensity'); 
    set(gca, 'FontSize', 14, 'TickDir','out','LineWidth',2);
    ax = gca; 
    xlim([12.5,14.3]);
    ylim([0,1]);
    %legend('PSNP blend', 'Fitted curve','Population 1','PSD population 1','Population 2','PSD population 2','Population 3','PSD population 3');
    %legend('PSNP blend', 'Fitted curve','Population 1','PSD population 1','Population 2','PSD population 2');
    %legend ('boxoff');
    %legend ('location','northwest');
    [fit_h, fit_loc, fit_FWHM] = findpeaks(y_fit,x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend','NPeaks',2);
    fit_para = [fit_h,fit_loc,fit_FWHM/2.35, polyval(p_var_2,fit_loc)]
    %findpeaks(yhat(:,1),x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend')
    [p1_h, p1_loc, p1_FWHM] = findpeaks(yhat(:,1),x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend');
    [p2_h, p2_loc, p2_FWHM] = findpeaks(yhat(:,2),x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend');
    [p3_h, p3_loc, p3_FWHM] = findpeaks(yhat(:,3),x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend');
    p1_para = [p1_h, p1_loc, p1_FWHM/2.35, polyval(p_var_2,p1_loc)]
    p2_para = [p2_h, p2_loc, p2_FWHM/2.35, polyval(p_var_2,p2_loc)]
    p3_para = [p3_h, p3_loc, p3_FWHM/2.35, polyval(p_var_2,p3_loc)]

elseif n == 4
    ft = ft4;
    yhat1 = (fa.a1 / (fa.c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b1) .^2) / (fa.c1 .^ 2));
    yhat2 = (fa.a2 / (fa.c2 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b2) .^2) / (fa.c2 .^ 2));
    yhat3 = (fa.a3 / (fa.c3 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b3) .^2) / (fa.c3 .^ 2));
    yhat4 = (fa.a4 / (fa.c4 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b4) .^2) / (fa.c4 .^ 2));
    yhat =[yhat1 yhat2 yhat3 yhat4];
    y_fit = sum(yhat,2);
    figure
    hold on
    plot(x,signal, 'LineWidth', 2,'color','#0D2747');
    plot(x,y_fit, 'LineWidth', 2,'color','#D06D91');       % fitted data
    plot(x,yhat(:,2), 'LineWidth', 2,'color','#567964','LineStyle','--');    % population 1
    plot(time,cor_chrom_1,'LineWidth',2,'color', '#567964'); % corrected chrom 1
    plot(x,yhat(:,3), 'LineWidth', 2,'color','#A7A7A7','LineStyle','--');    % population 2
    plot(time,cor_chrom_2,'LineWidth',2,'color', '#A7A7A7'); % corrected chrom 2
    plot(x,yhat(:,4), 'LineWidth', 2,'color','#9877AB','LineStyle','--');    % population 3
    plot(time,cor_chrom_3,'LineWidth',2,'color', '#9877AB'); % corrected chrom 3
    %plot(x,yhat(:,1), 'LineWidth', 2,'color','#CC9900','LineStyle','--');    % population 4
    hold off
    xlabel('Time (min)'); ylabel('Normalized intensity'); 
    set(gca, 'FontSize', 14, 'TickDir','out','LineWidth',2);
    ax = gca; 
    xlim([12.5,14.3]);
    ylim([0,1]);
    legend('PSNP blend', 'Fitted curve','Population 1','PSD population 1','Population 2','PSD population 2','Population 3','PSD population 3','Population 4');
    legend ('boxoff');
    legend ('location','northwest');
    [fit_h, fit_loc, fit_FWHM] = findpeaks(y_fit,x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend','NPeaks',2);
    fit_para = [fit_h,fit_loc,fit_FWHM/2.35, polyval(p_var_2,fit_loc)]
    [p1_h, p1_loc, p1_FWHM] = findpeaks(yhat(:,1),x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend');
    [p2_h, p2_loc, p2_FWHM] = findpeaks(yhat(:,2),x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend');
    [p3_h, p3_loc, p3_FWHM] = findpeaks(yhat(:,3),x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend');
    [p4_h, p4_loc, p4_FWHM] = findpeaks(yhat(:,4),x,'MinPeakHeight',0.1, 'Annotate', 'extents', 'WidthReference', 'halfheight','SortStr','descend');
    p1_para = [p1_h, p1_loc, p1_FWHM/2.35, polyval(p_var_2,p1_loc)]
    p2_para = [p2_h, p2_loc, p2_FWHM/2.35, polyval(p_var_2,p2_loc)]
    p3_para = [p3_h, p3_loc, p3_FWHM/2.35, polyval(p_var_2,p3_loc)]
    p4_para = [p4_h, p4_loc, p4_FWHM/2.35, polyval(p_var_2,p4_loc)]

elseif n == 5
    ft = ft5;
    yhat1 = (fa.a1 / (fa.c1 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b1) .^2) / (fa.c1 .^ 2));
    yhat2 = (fa.a2 / (fa.c2 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b2) .^2) / (fa.c2 .^ 2));
    yhat3 = (fa.a3 / (fa.c3 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b3) .^2) / (fa.c3 .^ 2));
    yhat4 = (fa.a4 / (fa.c4 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b4) .^2) / (fa.c4 .^ 2));
    yhat5 = (fa.a5 / (fa.c5 * sqrt(2 * pi()))) * exp(-0.5 * ((x - fa.b5) .^2) / (fa.c5 .^ 2));
    yhat =[yhat1 yhat2 yhat3 yhat4 yhat5];
    y_fit = sum(yhat,2);
    figure
    hold on
    plot(x,signal, 'LineWidth', 2,'color','#0D2747');
    plot(x,y_fit, 'LineWidth', 2,'color','#D06D91');       % fitted data
    plot(x,yhat(:,1), 'LineWidth', 2,'color','#567964','LineStyle','--');    % population 1
    plot(x,yhat(:,2), 'LineWidth', 2,'color','#A7A7A7','LineStyle','--');    % population 2
    plot(x,yhat(:,3), 'LineWidth', 2,'color','#9877AB','LineStyle','--');    % population 3
    plot(x,yhat(:,4), 'LineWidth', 2,'color','#CC9900','LineStyle','--');    % population 4
    plot(x,yhat(:,5), 'LineWidth', 2,'color','#2A60A5','LineStyle','--');    % population 5
    hold off
    xlabel('Time (min)'); ylabel('Normalized intensity'); 
    set(gca, 'FontSize', 14, 'TickDir','out','LineWidth',2); 
    ax = gca; 
    xlim([12.5,14.3]);
    %ylim([0,1.1]);
    legend('PSNP blend', 'Fitted curve','Population 1','Population 2','Population 3','Population 4', 'Population 5');
    legend ('boxoff');
    legend ('location','northwest');

else
    error("ERROR: too many peaks for fitting, please increase treshold")
end


y_fit = sum(yhat,2);


areacoverage = 100-((sum(signal)-sum(y_fit))/sum(signal))*100;

%p_var_2 = [-0.00352544207707122	0.0976618734670780	-0.668567622117514];




end