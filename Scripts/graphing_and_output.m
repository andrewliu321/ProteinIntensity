%% Graph all measured intensity

%define number of images analyzed
numofimages = size(Collect_APi_index_sum_Fat,2);
    
%choose smoothing parameter for moving line average
span=100;    



%% scatter graph 
numofimages = size(Collect_APi_index_sum_Fat,2);

% (5th column = intensity)
    f51=figure(51);
    clf
    f41=figure(41);
    clf
for pos = 1:numofimages %for now this will only plot one, but could plot a group of samples together
    f51=figure(51);
    sgtitle({'Erased Peripodial Sum Surface Projected','DV Boundary'},'FontSize',30)
    subplot(1,2,1)
    title('Fat GFP','FontSize',20)
    hold on
    scatter(Collect_APi_index_sum_Fat{pos}(:,4)/pixel_per_um,Collect_APi_index_sum_Fat{pos}(:,3),2) 
    xlabel('Position Along the DV Boundary (um)','FontSize',15)
    ylabel('Relative Expression of Fat','FontSize',15)

    f41=figure(41);
    sgtitle({'Erased Peripodial Sum Surface Projected','AP Boundary'},'FontSize',30)
    subplot(1,2,1)
    title('Fat GFP','FontSize',20)
    hold on
    scatter(Collect_DVi_index_sum_Fat{pos}(:,4)/pixel_per_um,Collect_DVi_index_sum_Fat{pos}(:,3),2)
    xlabel('Position Along the AP Boundary (um)','FontSize',15)
    ylabel('Relative Expression of Fat','FontSize',15)

    f51=figure(51);
    subplot(1,2,2)
    title('anti Ds','FontSize',20)
    hold on
    scatter(Collect_APi_index_sum_Ds{pos}(:,4)/pixel_per_um,Collect_APi_index_sum_Ds{pos}(:,3),2)
    xlabel('Position Along the DV Boundary (um)','FontSize',15)
    ylabel('Relative Expression of Dachsous','FontSize',15)

    f41=figure(41);
    subplot(1,2,2)
    title('anti Ds','FontSize',20)
    hold on
    scatter(Collect_DVi_index_sum_Ds{pos}(:,4)/pixel_per_um,Collect_DVi_index_sum_Ds{pos}(:,3),2)
    xlabel('Position Along the AP Boundary (um)','FontSize',15)
    ylabel('Relative Expression of Dachsous','FontSize',15)

   
end

%% smooth graph
% column 1 and 2 = x,y coord
% colume 3 = raw intensity 
% column 4 = distance (from imsane)
% column 5 = smooth intensity (moving line average)

    f52=figure(52);
    clf
    f42=figure(42);
    clf
    
for pos = 1:numofimages

    Collect_APi_index_sum_Fat{pos}(:,5)=smooth(Collect_APi_index_sum_Fat{pos}(:,3),span);
    Collect_DVi_index_sum_Fat{pos}(:,5)=smooth(Collect_DVi_index_sum_Fat{pos}(:,3),span);
    
    f52=figure(52);
    sgtitle({'Smoothed Erased Peripodial Sum Surface Projected','DV Boundary'},'FontSize',30)
    subplot(1,2,1)
    title('Fat GFP','FontSize',20)
    hold on
    plot(Collect_APi_index_sum_Fat{pos}(:,4)/pixel_per_um,Collect_APi_index_sum_Fat{pos}(:,5),'-r','LineWidth',2)
    xlabel('Position Along the DV Boundary (um)','FontSize',15)
    ylabel('Relative Expression of Fat','FontSize',15)

    f42=figure(42);
    sgtitle({'Smoothed Erased Peripodial Sum Surface Projected','AP Boundary'},'FontSize',30)
    subplot(1,2,1)
    title('Fat GFP','FontSize',20)
    hold on
    plot(Collect_DVi_index_sum_Fat{pos}(:,4)/pixel_per_um,Collect_DVi_index_sum_Fat{pos}(:,5),'-r','LineWidth',2) 
    xlabel('Position Along the AP Boundary (um)','FontSize',15)
    ylabel('Relative Expression of Fat','FontSize',15)

end

for pos = 1:numofimages

    Collect_APi_index_sum_Ds{pos}(:,5)=smooth(Collect_APi_index_sum_Ds{pos}(:,3),span);
    Collect_DVi_index_sum_Ds{pos}(:,5)=smooth(Collect_DVi_index_sum_Ds{pos}(:,3),span);

    f52=figure(52);
    subplot(1,2,2)
    title('anti-Ds','FontSize',20)
    hold on
    plot(Collect_APi_index_sum_Ds{pos}(:,4)/pixel_per_um,Collect_APi_index_sum_Ds{pos}(:,5),'-r','LineWidth',2)
    xlabel('Position Along the DV Boundary (um)','FontSize',15)
    ylabel('Relative Expression of Dachsous','FontSize',15)

    f42=figure(42);
    subplot(1,2,2)
    title('anti-Ds','FontSize',20)
    hold on
    plot(Collect_DVi_index_sum_Ds{pos}(:,4)/pixel_per_um,Collect_DVi_index_sum_Ds{pos}(:,5),'-r','LineWidth',2)
    xlabel('Position Along the AP Boundary (um)','FontSize',15)
    ylabel('Relative Expression of Dachsous','FontSize',15)

end


%% normalize absolute x distance, 
% 6th column = distance index with respect to center
% 7th column = distance index with respect to total distance 

    f53=figure(53);
    clf
    f43=figure(43);
    clf
for pos = 1:numofimages


    Collect_APi_index_sum_Fat{pos}(:,6)= normalize_x_with_center(Collect_APi_index_sum_Fat{pos}(:,4));
    Collect_DVi_index_sum_Fat{pos}(:,6)= normalize_x_with_center(Collect_DVi_index_sum_Fat{pos}(:,4));

    
    f53=figure(53);
    sgtitle({'Normalized Dist Smoothed Sum Surface Project Measurements','DV Boundary','Normalized with center'},'FontSize',30)
    subplot(1,2,1)
    title('Fat GFP','FontSize',20)
    hold on
    plot(Collect_APi_index_sum_Fat{pos}(:,6),Collect_APi_index_sum_Fat{pos}(:,5),'LineWidth',2) 
    xlabel('Position Along the DV Boundary','FontSize',15)
    ylabel('Relative Expression of Fat','FontSize',15)
    
    f43=figure(43);
    title({'Normalized Dist Smoothed Sum Surface Project Measurements','AP Boundary','Normalized with center'},'FontSize',30)
    subplot(1,2,1)
    title('Fat GFP','FontSize',20)
    hold on
    plot(Collect_DVi_index_sum_Fat{pos}(:,6),Collect_DVi_index_sum_Fat{pos}(:,5),'LineWidth',2)
    xlabel('Position Along the AP Boundary','FontSize',15)
    ylabel('Relative Expression of Fat','FontSize',15)

        

end

for pos = 1:numofimages

    Collect_APi_index_sum_Ds{pos}(:,6)= normalize_x_with_center(Collect_APi_index_sum_Ds{pos}(:,4));
    Collect_DVi_index_sum_Ds{pos}(:,6)= normalize_x_with_center(Collect_DVi_index_sum_Ds{pos}(:,4));
    
       
    f53=figure(53);
    sgtitle({'Normalized Dist Smoothed Sum Surface Project Measurements','DV Boundary','Normalized with center'},'FontSize',30)
    subplot(1,2,2)
    title('anti-Ds','FontSize',20)
    hold on
    plot(Collect_APi_index_sum_Ds{pos}(:,6),Collect_APi_index_sum_Ds{pos}(:,5),'LineWidth',2) 
    xlabel('Position Along the DV Boundary','FontSize',15)
    ylabel('Relative Expression of Dachsous','FontSize',15)
    
    f43=figure(43);
    title({'Normalized Dist Smoothed Sum Surface Project Measurements','AP Boundary','Normalized with center'},'FontSize',30)
    subplot(1,2,2)
    title('anti-Ds','FontSize',20)
    hold on
    plot(Collect_DVi_index_sum_Ds{pos}(:,6),Collect_DVi_index_sum_Ds{pos}(:,5),'LineWidth',2)
    xlabel('Position Along the AP Boundary','FontSize',15)
    ylabel('Relative Expression of Dachsous','FontSize',15)

end
    

    
%% save figures
if not(exist([folder,'analysis']))
    mkdir([folder,'analysis'])
end


savefig(f51,[folder,'\analysis\','Scatter DV'])
saveas(f51, [folder,'\analysis\','Scatter DV.jpg'])

savefig(f41,[folder,'\analysis\','Scatter AP'])
saveas(f41, [folder,'\analysis\','Scatter AP.jpg'])

savefig(f52,[folder,'\analysis\','Smoothed line DV'])
saveas(f52, [folder,'\analysis\','Smoothed line DV.jpg'])

savefig(f42,[folder,'\analysis\','Smoothed line AP'])
saveas(f42, [folder,'\analysis\','Smoothed line AP.jpg'])

f53 = figure(53);
savefig(f53,[folder,'\analysis\','Normalized Dist Smoothed Line DV'])
saveas(f53, [folder,'\analysis\','Normalized Dist Smoothed Line DV.jpg'])

f43 = figure(43);
savefig(f43,[folder,'\analysis\','Normalized Dist Smoothed Line AP'])
saveas(f43, [folder,'\analysis\','Normalized Dist Smoothed Line AP.jpg'])

%% save to csv

Fat_AP_Output = export_to_csv(Collect_APi_index_sum_Fat,pixel_per_um);
csvwrite([folder,'analysis\measured AP Fat in um.csv'],Fat_AP_Output) 

Fat_DV_Output = export_to_csv(Collect_DVi_index_sum_Fat,pixel_per_um);
csvwrite([folder,'analysis\measured DV Fat in um.csv'],Fat_DV_Output) 

Ds_AP_Output = export_to_csv(Collect_APi_index_sum_Ds,pixel_per_um);
csvwrite([folder,'analysis\measured AP Ds in um.csv'],Ds_AP_Output) 

Ds_DV_Output = export_to_csv(Collect_DVi_index_sum_Ds,pixel_per_um);
csvwrite([folder,'analysis\measured DV Ds in um.csv'],Ds_DV_Output) 
     
     
%% normalize X again:


normalized_Fat_AP_output = export_to_csv_normalized(Fat_AP_Output);
csvwrite([folder,'analysis\measured AP Fat normalized.csv'],normalized_Fat_AP_output) 

normalized_Fat_DV_output = export_to_csv_normalized(Fat_DV_Output);
csvwrite([folder,'analysis\measured DV Fat normalized.csv'],normalized_Fat_DV_output) 


normalized_Ds_AP_output = export_to_csv_normalized(Ds_AP_Output);
csvwrite([folder,'analysis\measured AP Ds normalized.csv'],normalized_Ds_AP_output) 

normalized_Ds_DV_output = export_to_csv_normalized(Ds_DV_Output);
csvwrite([folder,'analysis\measured DV Ds normalized.csv'],normalized_Ds_DV_output) 
          
          
          
          
          
%% functions:

function final_output = export_to_csv(Collect_index,pixel_per_um)
                   
        numofimages = size(Collect_index,2);

     for pos = 1:numofimages
        size_Ant(pos) = min(Collect_index{pos}(:,4));
        size_Pos(pos) = max(Collect_index{pos}(:,4));
         
     end
    
     smallest_ant = round(min(size_Ant)/pixel_per_um);
     biggest_pos = round(max(size_Pos)/pixel_per_um);
     biggest_size = (biggest_pos-smallest_ant);
     
     
     first_column = transpose(smallest_ant:biggest_pos);
     
     final_output = NaN(biggest_size+1,numofimages+1);
     
     final_output(:,1)=first_column;
     
     for pos = 1:numofimages
             
         x= Collect_index{pos}(:,4)/pixel_per_um;
         v = Collect_index{pos}(:,5);
         
         most_ant = round(Collect_index{pos}(1,4)/pixel_per_um);
         most_pos = round(Collect_index{pos}(end,4)/pixel_per_um);
         
         xq = (most_ant:most_pos);
         
         vq = interp1(x,v,xq);
         
         start_index = find(first_column==most_ant);
         end_index = find(first_column==most_pos);
                  
         final_output(start_index:end_index,pos+1) = vq;        
     end
end
          
function normalized_output = export_to_csv_normalized(A)

     %normalize
     for i = 2:size(A,2)
         
        TF = isnan(A(:,i));
        first_not_nan_index = find(TF==0,1);
        second_not_nan_index = find(TF==0,1,'last'); 
        
        x2 = A(first_not_nan_index,1):A(second_not_nan_index,1);
        zeroth_index = find(x2 == 0);

        neg_1_to_0 = x2(1:zeroth_index)./abs(x2(1));
        zero_to_pos_1 = x2(zeroth_index+1:end)./x2(end);

        neg_1_to_1 = horzcat(neg_1_to_0,zero_to_pos_1);
        
        x2q = -1:0.02:1;
        x2q=transpose(x2q);
        v2 = A(first_not_nan_index:second_not_nan_index,i);
        v2 = transpose(v2);
        
        v2q = interp1(neg_1_to_1,v2,x2q);
        
        normalized_output(:,1) = x2q;
        normalized_output(:,i)=v2q;

        figure(1)
        hold on
        plot(normalized_output(:,1),normalized_output(:,i))
         
     end
          
end


function Y = normalize_int_by_center(measured_intensity)
%find the center using 4th column (reg index)
center=find(measured_intensity(:,4)==0);
X=measured_intensity(1:center,5);

[n,m] = size(measured_intensity(:,5));
  Y = zeros(n,1);
  xmin = measured_intensity(center,5); %min(X);
  xmax = max(X);
  for jj = 1:center
    Y(jj) = (X(jj) - xmin) / (xmax-xmin);
  end
  
  X2=measured_intensity(center:end,5);
  xmax2 = max(X2);
  for jj = center:n
    Y(jj) = (measured_intensity(jj,5) - xmin) / (xmax2-xmin);
  end
  
end


function minus_min_int = minus_min(intensity)
% subtract by minimum for intensity so its all at 0

minimum = min(intensity);
minus_min_int = intensity-minimum;

end

function max_is_1 = divide_max1(intensity)
% divide by maximum so it turns to 1
% maximum defined by regular maximum

maximum = max(intensity);
max_is_1 = intensity./maximum;

end

function max_is_1 = divide_max2(index,intensity)
% divide by maximum so it turns to 1
% maximum defined by average of the to extreme indexes

zeroth_index = find(index == 0);

maximum1 = max(intensity(zeroth_index:end));
maximum2 = max(intensity(1:zeroth_index));
maximum = (maximum1+maximum2)/2;
max_is_1 = intensity./maximum;

end


function neg_1_to_1 = normalize_x_with_center(min_to_max)
% turn the index from -max to 0 to max
% to -1 to 0 to 1

zeroth_index = find(min_to_max == 0);

neg_1_to_0 = min_to_max(1:zeroth_index,1)./abs(min_to_max(1,1));
zero_to_pos_1 = min_to_max(zeroth_index+1:end,1)./min_to_max(end,1);

neg_1_to_1 = vertcat(neg_1_to_0,zero_to_pos_1);

end
    
    
function column2 = normalize_x_all(column1)
% turn the index from -max to 0 to max
% to -1 to 0 to 1

smallest=min(column1);
column1 = column1 + abs(smallest);
column2 = column1./(column1(end,1));

end
    
    
    
function [] = plotGFP(DataSetAP,DataSetDV, positions2plot, color, xlimAP, ylimAP, xlimDV, ylimDV,figNum, totalBinNum, BinNum, YonLR, titles)
% DataSetAP = which dataset to plot from, ie sorted_APi_index_sum
% DataSetDV = which dataset to plot from, ie sorted_DVi_index_sum
% positions2plot = which ones to plot ie [1:4,6]
% clor = which color lines
% xlimAP = x limit for AP axis ([-250 250])
% ylimAP = x limit for AP axis ([0 250])
% xlimDV = x limit for DV axis ([-120 80])
% ylimDV = x limit for DV axis ([0 250])
% figNum = which figure number ie 10
% totalBinNuum = how many bins we are plotting, 200-250, 250-300, 300-350,
% WPP would be 4
% BinNum, which one are we plotting now 200-250 = 1
% YonLR = are we plotting multiple genotypes? 0 = yes, 1 = Left,2 = right
% titles = {'200-250'}


    for pos = positions2plot
        figure(figNum)
        subplot(2,totalBinNum,BinNum)
        hold on
        if YonLR == 1
            yyaxis left
        elseif YonLR == 2
            yyaxis right
        end
        plot(DataSetAP{pos}(:,4)/pixel_per_um,DataSetAP{pos}(:,5),color,'LineWidth',2)%'-k''Color',color_list_wt(pos,:)
        ylim(ylimAP)
        xlim(xlimAP)
        title(titles,'FontSize',20)

        figure(figNum)
        subplot(2,totalBinNum,BinNum+totalBinNum)
        hold on
        if YonLR == 1
            yyaxis left
        elseif YonLR == 2
            yyaxis right
        end
        plot(DataSetDV{pos}(:,4)/pixel_per_um,DataSetDV{pos}(:,5),color,'LineWidth',2)%'-k''Color',color_list_wt(pos,:),
        ylim(ylimDV)
        xlim(xlimDV)
    end
end
