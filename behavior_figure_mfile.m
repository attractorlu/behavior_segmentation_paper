
%% plot 2-5 comparation both mice and monkey
clear
colors = [...
    255/255 224/255 170/255; ...
    134/255 125/255 176/255; ...
    108/255 162/255 153/255; ...
  %  0 0 0; ...
];

category = {'2_cross_last_day', '5_cross_first_day',...
                      '2_iso_last_day', '5_iso_first_day',...
                      '2_iso_tex_last_day', '5_iso_tex_first_day'};
%load('c_7_iso.mat');
%c_iso=correct_rate_total_sum;


all = zeros(12,6);


for i=1:length(category)
    files_to_load = find_files( ['*' category{i} '*']);
    n = length(files_to_load);
    cr = [];
    for j=1:n
        tmp = load(files_to_load{j});
        cr = cat(1, cr, tmp.correct_rate_total_sum);
    end
    all(:,i)=cr;
end

num=sqrt(12);
for i =1:size(all,2)
    sem(i)=std(all(:,i))/num;
end 

for i =1:size(all,2)
    mean_num(i)=mean(all(:,i));
end
%plot(cr_mean,'-o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:),'LineWidth', 2);

load ('D:\experiment\mouse training\summary\2_5_compare\2_5_compare_mice.mat')

gray=[0.8 0.8 0.8];
symbor_size=8;

figure ( 'position',  [ 326   329    684   420 ])
%plot([1:2], all(:,1:2),'.-','Color',[0.5 0.5 0.5] ,'MarkerSize',20);%colors(1,:)
plot( [1:2], all(:,1:2),'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2, 'MarkerSize',symbor_size);

hold on
plot([0.6 1.4], [mean_num(1) mean_num(1)], 'Color', colors(1,:), 'MarkerSize',20, 'LineWidth', 3);
plot([1.6 2.4], [mean_num(2) mean_num(2)] ,'Color', colors(1,:), 'MarkerSize',20, 'LineWidth', 3);
%errorbar(1,mean_num(1),sem(1),'.','Color',colors(1,:), 'MarkerSize',20, 'LineWidth', 1);
%errorbar(2,mean_num(2),sem(2),'.','Color',colors(1,:), 'MarkerSize',20, 'LineWidth', 1);


hold on
plot([3:4], all(:,3:4),'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2,'MarkerSize',symbor_size);%colors(2,:)
hold on
plot([2.6 3.4], [mean_num(3) mean_num(3)],'Color', colors(2,:), 'MarkerSize',20, 'LineWidth', 3);
plot([3.6 4.4], [mean_num(4) mean_num(4)],'Color', colors(2,:), 'MarkerSize',20, 'LineWidth', 3);
%errorbar(3,mean_num(3),sem(3),'.','Color',colors(2,:), 'MarkerSize',20, 'LineWidth', 1);
%errorbar(4,mean_num(4),sem(4),'.','Color',colors(2,:), 'MarkerSize',20, 'LineWidth', 1);

hold on
plot([5:6], all(:,5:6),'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2, 'MarkerSize',symbor_size);%colors(3,:)
hold on
plot([4.6 5.4], [mean_num(5) mean_num(5)],'Color', colors(3,:), 'MarkerSize',20, 'LineWidth', 3);
plot([5.6 6.4], [mean_num(6) mean_num(6)],'Color', colors(3,:), 'MarkerSize',20, 'LineWidth', 3);
%errorbar(5,mean_num(5),sem(5), '.','Color',colors(3,:), 'MarkerSize',20, 'LineWidth', 1);
%

load ('D:\experiment\ephys in mouse\paper_draft\Janis_Monkey_ObjectDetection\2_5_compare_monkey.mat')
% monkey data
gray_monkey=[0.5 0.5 0.5];
symbor_size_monkey=10;
hold on   % data  transposition
plot([8 9], [monkey_all(1,1)*100 monkey_all(1,2)*100],'-s', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:),'LineWidth', 2, 'MarkerSize',symbor_size_monkey);
plot([10 11], [monkey_all(1,3)*100 monkey_all(1,4)*100],'-s', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:),'LineWidth', 2, 'MarkerSize',symbor_size_monkey);
plot([12 13], [monkey_all(1,5)*100 monkey_all(1,6)*100],'-s', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:),'LineWidth', 2, 'MarkerSize',symbor_size_monkey);

rand=0.3;


plot([8+rand 9+rand], [monkey_all(2,1)*100 monkey_all(2,2)*100],'-^', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:),'LineWidth', 2, 'MarkerSize',symbor_size_monkey);
plot([10+rand 11+rand], [monkey_all(2,3)*100 monkey_all(2,4)*100],'-^', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:),'LineWidth', 2, 'MarkerSize',symbor_size_monkey);
plot([12+rand 13+rand], [monkey_all(2,5)*100 monkey_all(2,6)*100],'-^', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:),'LineWidth', 2, 'MarkerSize',symbor_size_monkey);
% 
% plot([9], monkey_all(:,2)'*100,'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2, 'MarkerSize',10);
% plot([10], monkey_all(:,3)'*100,'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2, 'MarkerSize',10);
% plot([11], monkey_all(:,4)'*100,'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2, 'MarkerSize',10);
% plot([12], monkey_all(:,5)'*100,'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2, 'MarkerSize',10);
% plot([13], monkey_all(:,6)'*100,'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2, 'MarkerSize',10);




ylim([40 100]);
xlim([0 13.5]);

xlabel('Days');
ylabel('Performance (% Correct)');
 set(gca,'XTick', 1:1:13 );
   
box off 
% legend('mouse 1', 'mouse 2', 'mouse 3', 'mouse 4');
set(gca, 'LineWidth', 2, 'FontSize', 14 ,'TickDir','out');
%legend('cross grating', 'iso grating', 'iso texture');

%% learning rate of cross iso and nat
clear
category = {'2_cross_grating', '2_iso_grating', '2_iso_texture'};
%category = {'2_cross_grating', '2_iso_grating', '2_iso_texture'};%,'g10_b'};
 %category = {'5_cross_grating', '5_iso_grating', '5_iso_texture'};
 %category = {'30_cross', '30_iso'};
 %category = {'2_iso_grating'};
 %category = {'5_cross_grating', '5_iso_grating', '5_iso_texture'};
 %category = {'30_cross', '30_iso'};
  %category = {'2_cross_texture', '2_iso_texture'};
  %category = {'g10_b'};
  
colors = [...
   % 0.8 0 0.8;
  % 0.3 0 0.3;...
     255/255 224/255 170/255; ...
    134/255 125/255 176/255; ...
    108/255 162/255 153/255; ...
   % 0 0 0; ...
];

colors_two = [...
   % 0.8 0 0.8;
  % 0.3 0 0.3;...
     244/255 242/255 189/255; ...
    198/255 189/255 240/255; ...
    172/255 226/255 217/255; ...
    %0.5 0.5 0.5; ...
];

figure ( 'position',  [ 326   329   684   420 ])
for i=1:length(category)
    files=find_files( ['*' category{i} '*']);
    n = length(files);
    for j=1:n
        data(j) = load(files{j});
        nday(j)= size(data(j).correct_rate_total_sum,2);
    end
    nday_min = min(nday);
    cr = [];
    for j=1:n
        cr = cat(1, cr, data(j).correct_rate_total_sum(:,1:nday_min));
    end
    cr_mean = mean(cr);
    cr_sem = std(cr)/sqrt(size(cr,1));
    %error_area(1:nday_min, cr_mean, cr_sem, min(colors(i,:)+0.5, [1 1 1]));
    error_area(1:nday_min, cr_mean, cr_sem, colors_two(i,:));
    hold on
    line(i) = plot(cr_mean,'-o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:),'LineWidth', 2, 'MarkerSize',10);
    
    ylim([40 100]);
    xlim([0.8 13+ 0.2]);
    
    %xlim([0.8 nday_min + 0.2]);
   %set(gca,'XTick', 1:1:nday_min );
    set(gca,'XTick', 1:1:14+0.2 );
   

xlabel('Days');
ylabel('Performance (% Correct)');
% legend('mouse 1', 'mouse 2', 'mouse 3', 'mouse 4');
set(gca, 'LineWidth', 2, 'FontSize', 14, 'TickDir','out');
    
end
legend(line, strrep(category, '_', ' '))
%%  plot static comparation
load ('D:\experiment\ephys in mouse\paper_draft\Janis_Monkey_ObjectDetection\two_monkey.mat')

oz_cross_last_trial = oz.LearningCurveOnSquareCrossV1(end);
oz_cross_first_static = oz.LearningCurveOnSquareCrossStatic(1);

oz_iso_last_trial = oz.LearningCurveOnSquareIsoV1(end);
oz_iso_first_static = oz.LearningCurveOnSquareIsoStatic(1);

oz_nat_last_trial = oz.LearningCurveOnSquareNatV1(end);
oz_nat_first_static = oz.LearningCurveOnSquareNatStatic(1);


Al_cross_last_trial = Al.LearningCurveOnSquareCrossV1(end);
Al_cross_first_static = Al.LearningCurveOnSquareCrossStatic(1);

Al_iso_last_trial = Al.LearningCurveOnSquareIsoV1(end);
Al_iso_first_static = Al.LearningCurveOnSquareIsoStatic(1);

Al_nat_last_trial = Al.LearningCurveOnSquareNatV1(end);
Al_nat_first_static = Al.LearningCurveOnSquareNatStatic(1);

monkey_static(1,: ) = [oz_cross_last_trial    Al_cross_last_trial   oz_iso_last_trial   Al_iso_last_trial   oz_nat_last_trial   Al_nat_last_trial];
monkey_static(2,: ) = [oz_cross_first_static   Al_cross_first_static   oz_iso_first_static  Al_iso_first_static  oz_nat_first_static  Al_nat_first_static];





c_last_cross_grating = [87.395 84.27 83.408 89.666];
c_first_cross_grating_static = [84.524 85.065 85.256 90];

e_last_cross_grating = [82.985 82.985 92.667 92.308];  
e_first_cross_grating_static = [87.671 81.395 88.571 89.524];

c_last_iso_grating = [79.636 64.807 63.462 79.412];
c_first_iso_grating_static = [89.172 69.369 65.812 79.31];

e_last_iso_grating =[67.978 59.136 69.157 71.314];  
e_first_iso_grating_static = [79.231 69.582 81.271 86.392];


b_last_nat = [91.398 91.437 65.882 83.991];  
b_first_nat_static= [88.83  89.23 57.14 82.54];

e_last_nat = [84.74 84.57 79.11 72.67];
e_first_nat_static = [73.308 86.99 83.668 72.274]; 




% mouse_static(:,1 ) = [c_last_cross_grating'];% e_last_cross_grating'];
% mouse_static(:,2 ) = [c_first_cross_grating_static'];% e_first_cross_grating_static'];
% mouse_static(:,3 ) = [c_last_iso_grating'];% e_last_iso_grating'];
% mouse_static(:,4 ) = [c_first_iso_grating_static'];% e_first_iso_grating_static'];
% 
% mouse_static(:,5 ) = [b_last_nat'];%;0;0;0;0 ];
% mouse_static(:,6 ) = [b_first_nat_static'];%;0;0;0;0 ];

%
mouse_static(:,1 ) = [e_last_cross_grating'; c_last_cross_grating' ];% e_last_cross_grating'];
mouse_static(:,2 ) = [e_first_cross_grating_static'; c_first_cross_grating_static'];% e_first_cross_grating_static'];
mouse_static(:,3 ) = [e_last_iso_grating'; c_last_iso_grating'];% e_last_iso_grating'];
mouse_static(:,4 ) = [e_first_iso_grating_static'; c_first_iso_grating_static'];% e_first_iso_grating_static'];

mouse_static(:,5 ) = [e_last_nat'; b_last_nat' ];%;0;0;0;0 ];
mouse_static(:,6 ) = [e_first_nat_static'; b_first_nat_static'];%;0;0;0;0 ];
%
  
colors = [...
   % 0.8 0 0.8;
  % 0.3 0 0.3;...
     255/255 224/255 170/255; ...
    134/255 125/255 176/255; ...
    108/255 162/255 153/255; ...
   % 0 0 0; ...
];


for i =1:6
    mean_num(i)=mean(mouse_static(:,i));
end

%mean_num(5) = mean(mouse_static(1:4, 5));
%mean_num(6) = mean(mouse_static(1:4, 6));


gray=[0.8 0.8 0.8];
symbor_size=8;

figure ( 'position',  [ 326   329   684   420 ])
%plot([1:2], all(:,1:2),'.-','Color',[0.5 0.5 0.5] ,'MarkerSize',20);%colors(1,:)
plot( [1:2], mouse_static(:,1:2),'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2, 'MarkerSize',symbor_size);

hold on
plot([0.6 1.4], [mean_num(1) mean_num(1)], 'Color', colors(1,:), 'MarkerSize',20, 'LineWidth', 3);
plot([1.6 2.4], [mean_num(2) mean_num(2)] ,'Color', colors(1,:), 'MarkerSize',20, 'LineWidth', 3);
%errorbar(1,mean_num(1),sem(1),'.','Color',colors(1,:), 'MarkerSize',20, 'LineWidth', 1);
%errorbar(2,mean_num(2),sem(2),'.','Color',colors(1,:), 'MarkerSize',20, 'LineWidth', 1);


hold on
plot([3:4], mouse_static(:,3:4),'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2, 'MarkerSize',symbor_size);%colors(2,:)
hold on
plot([2.6 3.4], [mean_num(3) mean_num(3)],'Color', colors(2,:), 'MarkerSize',20, 'LineWidth', 3);
plot([3.6 4.4], [mean_num(4) mean_num(4)],'Color', colors(2,:), 'MarkerSize',20, 'LineWidth', 3);
%errorbar(3,mean_num(3),sem(3),'.','Color',colors(2,:), 'MarkerSize',20, 'LineWidth', 1);
%errorbar(4,mean_num(4),sem(4),'.','Color',colors(2,:), 'MarkerSize',20, 'LineWidth', 1);

hold on
plot([5:6], mouse_static(:,5:6),'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2, 'MarkerSize',symbor_size);%colors(3,:)
hold on
plot([4.6 5.4], [mean_num(5) mean_num(5)],'Color', colors(3,:), 'MarkerSize',20, 'LineWidth', 3);
plot([5.6 6.4], [mean_num(6) mean_num(6)],'Color', colors(3,:), 'MarkerSize',20, 'LineWidth', 3);
%errorbar(5,mean_num(5),sem(5), '.','Color',colors(3,:), 'MarkerSize',20, 'LineWidth', 1);
%


gray_monkey=[0.5 0.5 0.5];
symbor_size_monkey=10;
hold on   % data  transposition

plot([8 9], [monkey_static(1,1)*100 monkey_static(2,1)*100],'-s', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:),'LineWidth', 2, 'MarkerSize',symbor_size_monkey);
plot([10 11], [monkey_static(1,3)*100 monkey_static(2,3)*100],'-s', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:),'LineWidth', 2, 'MarkerSize',symbor_size_monkey);
plot([12 13], [monkey_static(1,5)*100 monkey_static(2,5)*100],'-s', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:),'LineWidth', 2, 'MarkerSize',symbor_size_monkey);

rand=0.3;


plot([8+rand 9+rand], [monkey_static(1,2)*100 monkey_static(2,2)*100],'-^', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:),'LineWidth', 2, 'MarkerSize',symbor_size_monkey);
plot([10+rand 11+rand], [monkey_static(1,4)*100 monkey_static(2,4)*100],'-^', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:),'LineWidth', 2, 'MarkerSize',symbor_size_monkey);
plot([12+rand 13+rand], [monkey_static(1,6)*100 monkey_static(2,6)*100],'-^', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:),'LineWidth', 2, 'MarkerSize',symbor_size_monkey);
% plot monkey data



% 
% hold on   % data  transposition
% plot([8:9], monkey_static(:,1:2)*100,'-o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5],'LineWidth', 2);
% plot([10:11], monkey_static(:,3:4)*100,'-o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5],'LineWidth', 2);
% plot([12:13], monkey_static(:,5:6)*100,'-o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5],'LineWidth', 2);

ylim([30 100]);
xlim([0 14]);
xlabel('Days');
ylabel('Performance (% Correct)');
set(gca,'XTick', 1:1:13 );
box off 
% legend('mouse 1', 'mouse 2', 'mouse 3', 'mouse 4');
set(gca, 'LineWidth', 2, 'FontSize', 14 ,'TickDir','out');

%% plot the static differences

colors = [...
   % 0.8 0 0.8;
  % 0.3 0 0.3;...
     255/255 224/255 170/255; ...
    134/255 125/255 176/255; ...
    108/255 162/255 153/255; ...
   % 0 0 0; ...
];


diff_monkey=(monkey_static(2,:)-monkey_static(1, :));%./monkey_static(1,:);
diff_monkey = reshape(diff_monkey,2,3);

diff_mouse(:,1)= (mouse_static(:,2)-mouse_static(:,1));%./mouse_static(:,1);
diff_mouse(:,2)= (mouse_static(:,4)-mouse_static(:,3));%./mouse_static(:,3);
diff_mouse(:,3)= (mouse_static(:,6)-mouse_static(:,5));%./mouse_static(:,5);
%%
figure ( 'position',  [ 326   329  550   420 ])

for i=1:3
%plot(i, diff_mouse(:,i),'o', 'Color', [1 1 1], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize',15, 'LineWidth', 2);%colors(3,:)
scatter([i i i i i i i i], diff_mouse(:,i), 'filled','SizeData',150,'MarkerFaceColor', [0.5 0.5 0.5]);%colors(3,:)
hold on
end
alpha(0.3)
mean_mouse=mean(diff_mouse);
for i=1:3
plot([i-0.3 i+0.3],[ mean_mouse(:,i) mean_mouse(:,i)],'-', 'Color', colors(i,:), 'MarkerSize',20, 'LineWidth', 3);%colors(3,:)
hold on
end

rand=0.2;
hold on  % first row is oz, second is al
scatter([5 6 7], diff_monkey(1,:)*100, 's','filled','SizeData',200,'MarkerFaceColor', [0.5 0.5 0.5]);
scatter([5+rand 6+rand 7+rand], diff_monkey(2,:)*100, '^','filled','SizeData',200,'MarkerFaceColor', [0.5 0.5 0.5]);
% plot(5, diff_monkey(:,1),'o', 'Color', [1 1 1], 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',15,'LineWidth', 2);%colors(3,:)
% plot(6, diff_monkey(:,2),'o', 'Color', [1 1 1], 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',15,'LineWidth', 2);%colors(3,:)
% plot(7, diff_monkey(:,3),'o', 'Color', [1 1 1], 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',15,'LineWidth', 2);%colors(3,:)
alpha(0.5)
plot([0.5 7.5],[0 0],'--', 'LineWidth', 3, 'Color',[0 0 0]);

ylim([-70 20]);
xlim([0.5 7.5]);
xlabel('Days');
%ylabel('% change in performance when removing motion cues');
box off 
% legend('mouse 1', 'mouse 2', 'mouse 3', 'mouse 4');
set(gca, 'LineWidth', 2, 'FontSize', 14 ,'TickDir','out');
%%  plot performance of cross iso and nat in monkey
load ('D:\experiment\ephys in mouse\paper_draft\Janis_Monkey_ObjectDetection\Analysis Object Detection (1)\two_monkey_raw_data.mat')


colors = [...
   % 0.8 0 0.8;
  % 0.3 0 0.3;...
     255/255 224/255 170/255; ...
    134/255 125/255 176/255; ...
    108/255 162/255 153/255; ...
   % 0 0 0; ...
];



al_cross=sum( al.LearningCurveOnSquareCrossV1)/length(al.LearningCurveOnSquareCrossV1)*100;
oz_cross=sum( oz.LearningCurveOnSquareCrossV1)/length(oz.LearningCurveOnSquareCrossV1)*100;

al_iso=sum( al.LearningCurveOnSquareIsoV1)/length(al.LearningCurveOnSquareIsoV1)*100;
oz_iso=sum( oz.LearningCurveOnSquareIsoV1)/length(oz.LearningCurveOnSquareIsoV1)*100;

al_nat=sum( al.LearningCurveOnSquareNatV1)/length(al.LearningCurveOnSquareNatV1)*100;
oz_nat=sum( oz.LearningCurveOnSquareNatV1)/length(oz.LearningCurveOnSquareNatV1)*100;

figure ( 'position',  [ 326   329   421   420 ])
hold on
plot(0.9 , [al_cross],'^', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:),'MarkerSize',10, 'LineWidth', 2);
plot( 1.1, [oz_cross],'s', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:),'MarkerSize',10, 'LineWidth', 2);%colors(3,:)
plot([1.9 ], [al_iso],'^', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:),'MarkerSize',10, 'LineWidth', 2);
plot([2.1], [oz_iso],'s', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:),'MarkerSize',10, 'LineWidth', 2);%colors(3,:)
plot([2.9 ], [al_nat ],'^', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:),'MarkerSize',10, 'LineWidth', 2);
plot([3.1], [oz_nat],'s', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:),'MarkerSize',10, 'LineWidth', 2);%colors(3,:)

ylim([40 100]);
xlim([0.5 3.5]);
xlabel('Days');
%ylabel('% change in performance when removing motion cues');
box off 
% legend('mouse 1', 'mouse 2', 'mouse 3', 'mouse 4');
set(gca, 'LineWidth', 2, 'FontSize', 14 ,'TickDir','out');
%

%% learning curve of bgstatic and nat
clear

category = {'iso_texture', 'cross_texture_bgstatic'};  % cage d and e in plot
  %c iso texture; d_2 iso texture; e 2 iso texture;
  % c 2 cross texture bgstatic; d 2 cross texture bgstatic; e 2 cross
  % texture bgstatic;


colors = [...
   % 0.8 0 0.8;
  % 0.3 0 0.3;... 
  108/255 162/255 153/255; ...
     128/255 128/255 128/255; ...
   
   
   % 0 0 0; ...
];

colors_two = [...
   % 0.8 0 0.8;
  % 0.3 0 0.3;...    
  172/255 226/255 217/255; ...
      200/255 200/255 200/255;...
   

   % 0 0 0; ...
];

figure ( 'position',  [ 326   329   684   420 ])
for i=1:length(category)
    files=find_files( ['*' category{i} '*']);
    n = length(files);
    for j=1:n
        data(j) = load(files{j});
        nday(j)= size(data(j).correct_rate_total_sum,2);
    end
    nday_min = min(nday);
    cr = [];
    for j=1:n
        cr = cat(1, cr, data(j).correct_rate_total_sum(:,1:nday_min));
    end
    cr_mean = mean(cr);
    cr_sem = std(cr)/sqrt(size(cr,1));
    %error_area(1:nday_min, cr_mean, cr_sem, min(colors(i,:)+0.5, [1 1 1]));
    error_area(1:nday_min, cr_mean, cr_sem, colors_two(i,:));
    hold on
    line(i) = plot(cr_mean,'-o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:),'LineWidth', 2, 'MarkerSize',10);
    
    ylim([40 100]);
    xlim([0.8 13+ 0.2]);
    
    %xlim([0.8 nday_min + 0.2]);
   %set(gca,'XTick', 1:1:nday_min );
    set(gca,'XTick', 1:1:14+0.2 );
   

xlabel('Days');
ylabel('Performance (% Correct)');
% legend('mouse 1', 'mouse 2', 'mouse 3', 'mouse 4');
set(gca, 'LineWidth', 2, 'FontSize', 14, 'TickDir','out');
    
end
legend(line, strrep(category, '_', ' '))
%%
load('D:\experiment\mouse training\summary\7_test_texture\test_7.mat');
%cross b cross a iso b iso a
colors = [...
   % 0.8 0 0.8;
  % 0.3 0 0.3;...
     255/255 224/255 170/255; ...
    134/255 125/255 176/255; ...
    108/255 162/255 153/255; ...
   % 0 0 0; ...
];

mean_num=mean(test_7,1);
a=ones(12,4);
gray=[0.5 0.5 0.5];



figure ( 'position',  [ 326   329   684   420 ])
%plot([1:2], all(:,1:2),'.-','Color',[0.5 0.5 0.5] ,'MarkerSize',20);%colors(1,:)
plot( [1:2], test_7(:,1:2),'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2, 'MarkerSize',10);
%plot( [1:2], test_7(:,1:2),'-', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5],'LineWidth', 2, 'MarkerSize',10);
hold on

%scatter( [a(:,1) a(:,2)*2], test_7(:,1:2),'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5],'LineWidth', 2, 'MarkerSize',10)



hold on
plot([0.75 1.25], [mean_num(1) mean_num(1)], 'Color', [0 0 0], 'MarkerSize',20, 'LineWidth', 3);
plot([1.75 2.25], [mean_num(2) mean_num(2)] ,'Color',[0 0 0],  'MarkerSize',20, 'LineWidth', 3);
%errorbar(1,mean_num(1),sem(1),'.','Color',colors(1,:), 'MarkerSize',20, 'LineWidth', 1);
%errorbar(2,mean_num(2),sem(2),'.','Color',colors(1,:), 'MarkerSize',20, 'LineWidth', 1);

plot( [3:4], test_7(:,3:4),'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2, 'MarkerSize',10);

hold on
plot([2.75 3.25], [mean_num(3) mean_num(3)], 'Color', [0 0 0], 'MarkerSize',20, 'LineWidth', 3);
plot([3.75 4.25], [mean_num(4) mean_num(4)] ,'Color',[0 0 0], 'MarkerSize',20, 'LineWidth', 3);



%alpha(h, 0.5)
ylim([40 80]);
xlim([0 5]);
xlabel('Days');
ylabel('Performance (% Correct)');
box off 
% legend('mouse 1', 'mouse 2', 'mouse 3', 'mouse 4');
set(gca, 'LineWidth', 2, 'FontSize', 14 ,'TickDir','out');



%%  plot generalization of cross iso nat in mouse and monkey
load ('D:\experiment\ephys in mouse\paper_draft\Janis_Monkey_ObjectDetection\Analysis Object Detection (1)\two_monkey_raw_data.mat')

colors = [...
   % 0.8 0 0.8;
  % 0.3 0 0.3;...
     255/255 224/255 170/255; ...
    134/255 125/255 176/255; ...
    108/255 162/255 153/255; ...
   % 0 0 0; ...
];

al_cross=sum( al.LearningCurveOnSquareCrossGeneralized)/length(al.LearningCurveOnSquareCrossGeneralized)*100;
oz_cross=sum( oz.LearningCurveOnSquareCrossGeneralized)/length(oz.LearningCurveOnSquareCrossGeneralized)*100;

al_iso=sum( al.LearningCurveOnSquareIsoGeneralized)/length(al.LearningCurveOnSquareIsoGeneralized)*100;
oz_iso=sum( oz.LearningCurveOnSquareIsoGeneralized)/length(oz.LearningCurveOnSquareIsoGeneralized)*100;

al_nat=sum( al.LearningCurveOnSquareNatGeneralized)/length(al.LearningCurveOnSquareNatGeneralized)*100;
oz_nat=sum( oz.LearningCurveOnSquareNatGeneralized)/length(oz.LearningCurveOnSquareNatGeneralized)*100;

num=sqrt(2);
all_cross=[al_cross;oz_cross];
all_iso=[al_iso;oz_iso];
all_nat=[al_nat;oz_nat];
all=[all_cross all_iso all_nat ];

for i =1:3
    sem_monkey(i)=std(all(:,i))/num;
end 


monkey_cross_iso_nat=[mean([al_cross ;oz_cross]); mean([al_iso ;oz_iso]);mean([al_nat ;oz_nat])];

mouse_cross_iso_nat = [84.4974; 70.579; 54.8357];
sem_mouse=[1.7516;1.3737; 1.2921];



rand=[-0.05 0 0.05];
%%
figure ( 'position',  [ 326   329   684   420 ])
hold on

for i=1:3
   
  plot(1.3+rand(i), monkey_cross_iso_nat(i),'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:),'MarkerSize',8,'LineWidth', 2);  
  errorbar(1.3+rand(i),monkey_cross_iso_nat(i),sem_monkey(i),  'Color', colors(i,:),'LineWidth', 2)
  
end


for i=1:3
  plot(1, mouse_cross_iso_nat(i),'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:),'MarkerSize',8,'LineWidth', 2);  
  errorbar(1,mouse_cross_iso_nat(i),sem_mouse(i),  'Color', colors(i,:),'LineWidth', 2)
 
end
%errorbar(1,mean_num(1),sem(1),'.','Color',colors(1,:), 'MarkerSize',20, 'LineWidth', 1);

ylim([50 100]);
xlim([0.8 1.5]);
xlabel('Days');
ylabel('Performance (% Correct)');
box off 
% legend('mouse 1', 'mouse 2', 'mouse 3', 'mouse 4');
set(gca, 'LineWidth', 2, 'FontSize', 14 ,'TickDir','out');
%% plot luminance

mouse_lum_cagee = [92.283; 88.308; 90.338; 91.337];%cage e
mouse_lum_caged = [94.773; 85.714; 89.881; 83.468];
mouse_lum_cageb = [85.106; 76.712; 87.778; 90.164];% stage2_b_black_9.22.17
mouse_lum_cagec = [89.474; 95.522; 91.935; 85.246]; % c_black_9.25.17

mouse_lum = [mouse_lum_cagee;mouse_lum_caged; mouse_lum_cageb; mouse_lum_cagec ];

sem=std(mouse_lum)/sqrt(length(mouse_lum));



al_lum=sum( al.LearningCurveOnSquareLuminance)/length(al.LearningCurveOnSquareLuminance)*100;
oz_lum=sum( oz.LearningCurveOnSquareLuminance)/length(oz.LearningCurveOnSquareLuminance)*100;

monkey_lum=[al_lum; oz_lum];

figure ( 'position',  [ 326   329   684   420 ])
hold on
gray=[0.5 0.5 0.5];

   
 %plot([1 1 1 1 ], mouse_lum,'o', 'Color', gray, 'MarkerFaceColor', gray,'MarkerSize',8,'LineWidth', 2);  
 plot( 1, mean( mouse_lum),'o', 'Color', gray, 'MarkerFaceColor', gray,'MarkerSize',6,'LineWidth', 2); 
 errorbar(1, mean( mouse_lum),sem, 'Color', gray,'LineWidth', 2);
hold on
plot([1.3 ], monkey_lum(1),'^', 'Color', gray, 'MarkerFaceColor', gray,'MarkerSize',8,'LineWidth', 2);  
 plot([1.3 ], monkey_lum(2),'s', 'Color', gray, 'MarkerFaceColor', gray,'MarkerSize',8,'LineWidth', 2);  

ylim([50 100]);
xlim([0.5 1.5]);
xlabel('Days');
ylabel('Performance (% Correct)');
box off 
% legend('mouse 1', 'mouse 2', 'mouse 3', 'mouse 4');
set(gca, 'LineWidth', 2, 'FontSize', 14 ,'TickDir','out');


%% plot grating 2 noise of cage b
clear
files = find_files('*_b.csv');

files = files([2:end 1]);

correct_rate_all = [];
days = [];
days_total = 0;
for i = 1:length(files)
    
    trial_data = load_trial_data( files{i} );
    correct_rate = get_correct_rate( trial_data );
    correct_rate_all = [correct_rate_all correct_rate nan(4,1)];
    days = [ days   days_total+1:days_total+size(correct_rate, 2)  nan];
    days_total = days_total + size(correct_rate, 2);
%     ndays(i) = size(correct_rate, 2);
    
end

colors = [...
    253/255 51/255 51/255; ...
    41/255 231/255 91/255; ...
    41/255 141/255 255/255; ...
    0 0 0; ...
];


%figure('Position',[ 4         333        1913         420 ]); 

figure('Position',[ 294         471        1138         334 ]); 

hold on
for i = 1:4
    plot(days, correct_rate_all(i,:),'-o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:),'LineWidth', 1.5,'MarkerSize',7);
end


plot([0.5 33],[85 85],'--', 'Color', 'k', 'LineWidth', 2)

ylim([30 100]);
xlim([0.8 33]);
% set(gca,'XTick',[1:1:10]);

%xtick('TickDir', out)
xlabel('Days');
ylabel('Performance (% Correct)');
% legend('mouse 1', 'mouse 2', 'mouse 3', 'mouse 4');
set(gca, 'LineWidth', 2, 'FontSize', 14, 'TickDir','out');

%last day of stage10
% correct rate: 91.398; 91.437; 83.991
% num of trials: 558; 654; 431


new_texture=[58.97; 62.83; 63.72; 56.88]; 
%trial num: 234; 296;218 :1 2 4
% b_7_no_border_textures 
%D:\experiment\mouse training\data\fig_organized\b_7_textures\b_7_no_border_textures.2.15.csv

static_texture=[88.83; 89.23; 57.14; 82.54];
%trial num: 394 455 315
% D:\experiment\mouse training\data\fig_organized\no_gap_static_noise\g10_b_static.csv

hold on
for i=1:4
plot(29, new_texture(i),'-o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:),'MarkerSize',7)
end

hold on
for i=1:4
 plot(32, static_texture(i),'-o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:),'MarkerSize',7)   
end
%% compute significant of b gra2noise with new and static texture


%last day of stage10
% correct rate: 91.398; 91.437; 83.991
% num of trials: 558; 654; 431


new_texture=[58.97; 62.83; 63.72; 56.88]; 
%trial num: 234; 296;218 :1 2 4
% b_7_no_border_textures 
%D:\experiment\mouse training\data\fig_organized\b_7_textures\b_7_no_border_textures.2.15.csv

static_texture=[88.83; 89.23; 57.14; 82.54];
%trial num: 394 455 315  : 1 2 4

%
n1 =91;      % correct rate
N1 =654;  % trial numbers
n2 =89;       % correct rate
N2 =455;  % trial numbers

x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [ones(n1,1); repmat(2,N1-n1,1); ones(n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2);






%%  plot learning curve of cross iso and nat with raw data 

%  mean  response
clear 
 category = {'2_cross_grating', '2_iso_grating', '2_iso_texture'};
 %category = {'5_cross_grating', '5_iso_grating', '5_iso_texture'};
% category = {'30_cross', '30_iso'};
colors = [...
    %0.8 0 0.8;
   %0.3 0 0.3;...
    255/255 224/255 170/255; ...
    134/255 125/255 176/255; ...
    108/255 162/255 153/255; ...
  %  0 0 0; ...
];

figure ( 'position',  [ 326   329   684   420 ])
for i=1:length(category)
    files=find_files( ['*' category{i} '*']);
    n = length(files);
    for j=1:n
        data(j) = load(files{j});
        nday(j)= size(data(j).correct_rate_total_sum,2);
        for k=1:size(data(j).correct_rate_total_sum,1)
        plot(data(j).correct_rate_total_sum(k,:),'Color', colors(i,:)); 
        end% plot raw traces
        hold on
    end
    nday_min = min(nday);
    
%     for m=1:n
%         data(m) = load(files{m});
%         for k=1:size(data(m).correct_rate_total_sum,1)
%         plot(data(m).correct_rate_total_sum(k,nday_min)); 
%         end% plot raw traces
%         hold on
%     end
    
 
    
    cr = [];
    for j=1:n
        cr = cat(1, cr, data(j).correct_rate_total_sum(:,1:nday_min));
    end
    cr_mean = mean(cr);
    cr_sem = std(cr)/sqrt(size(cr,1));
    %error_area(1:nday_min, cr_mean, cr_sem, min(colors(i,:)+0.25, [1 1 1]));
    hold on
    line(i) = plot(cr_mean,'-o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:),'LineWidth', 2);
    
    ylim([40 100]);
    xlim([0.8 22+ 0.2]);
    
   % xlim([0.8 nday_min + 0.2]);
   set(gca,'XTick', 1:1:22 );
    %set(gca,'XTick', 1:1:14+0.2 );
   

xlabel('Days');
ylabel('Performance (% Correct)');
% legend('mouse 1', 'mouse 2', 'mouse 3', 'mouse 4');
set(gca, 'LineWidth', 2, 'FontSize', 14, 'TickDir','out');

box off

end

%% cumulative performance across trials in one session
colors = [...
    255/255 224/255 170/255; ...
    134/255 125/255 176/255; ...
    108/255 162/255 153/255; ...
  %  0 0 0; ...
];

cross_al = al.LearningCurveOnSquareCrossV1;
iso_al = al.LearningCurveOnSquareIsoV1;
nat_al = al.LearningCurveOnSquareNatV1;

cross_oz = oz.LearningCurveOnSquareCrossV1;
iso_oz = oz.LearningCurveOnSquareIsoV1;
nat_oz = oz.LearningCurveOnSquareNatV1;

% cross_al = al.LearningCurveOnSquareCrossGeneralized;
% iso_al = al.LearningCurveOnSquareIsoGeneralized;
% nat_al = al.LearningCurveOnSquareNatGeneralized;
% 
% cross_oz = oz.LearningCurveOnSquareCrossGeneralized;
% iso_oz = oz.LearningCurveOnSquareIsoGeneralized;
% nat_oz = oz.LearningCurveOnSquareNatGeneralized;

% compute correct rate

al_raw={cross_al;iso_al;nat_al};
oz_raw={cross_oz;iso_oz;nat_oz};

raw={al_raw; oz_raw};


figure ( 'position',  [ 326   329   684   420 ])
sym_type={'-','-'};

for m=1:2
for k =1:3
    n_trials= length(raw{m}{k});
    cum_cr = zeros(1, n_trials);
    temp_data=raw{m}{k};
for i = 1:n_trials
    temp = temp_data(1:i); 
    cum_cr(i)=sum(temp)/length(temp); 
    num = length(cum_cr);
    bin_num=10;
    cr_window=zeros(1,num);
    
    for j=bin_num:num
           cr_window(j)=mean(cum_cr(j-(bin_num-1):j));       
    end
    cr_window(1:bin_num)=cum_cr(1:bin_num);
    
  %  plot(cr_window*100,'-', 'Color', colors(k,:),'LineWidth', 4);
end
 plot(cr_window*100, 'LineStyle', sym_type{m}, 'Color', colors(k,:),'LineWidth', 4);
hold on
end
hold on
end

ylim([40 100]);
xlabel('Trials');
ylabel('Performance (% Correct)');
set(gca, 'LineWidth', 2, 'FontSize', 14, 'TickDir','out');
box off
%% 7 test cross and iso
colors = [...
    %0.8 0 0.8;
   %0.3 0 0.3;...
   0.5 0.5 0.5
    255/255 224/255 170/255; ...
    134/255 125/255 176/255; ...
    108/255 162/255 153/255; ...
  %  0 0 0; ...
];
gray=[0.5 0.5 0.5];


% cross_7_test = xlsread('D:\experiment\ephys in mouse\paper_draft\behavior_Figure\30texture\7\cross_7_test.xlsx'); % before after cage c d e
% iso_7_test = xlsread('D:\experiment\ephys in mouse\paper_draft\behavior_Figure\30texture\7\iso_7_test.xlsx');
% %delete d 3 4 
% cross_7_test=cross_7_test([1:4 7:end], :);  % c1-4, d 1-2, e 1-4
% iso_7_test=iso_7_test([1:4 7:end], :);  % c1-4, d 1-2, e 1-4

%delete cage c an cage d 3 4
cross_7_test = xlsread('D:\experiment\ephys in mouse\paper_draft\behavior_Figure\30texture\7\cross_7_test_delete_cage_c.xlsx'); % before after cage c d e
iso_7_test = xlsread('D:\experiment\ephys in mouse\paper_draft\behavior_Figure\30texture\7\iso_7_test_delete_cage_c.xlsx');

x = ones(6,1);

%

% x=ones(10,1); % c1-4, d 1-2, e 1-4


mean_cross=mean(cross_7_test);
mean_iso=mean(iso_7_test);


% plot cross baseline training test
%%
figure ( 'position',  [ 326         329        1130         420 ])
scatter( x*(-2), cross_7_test(:,1),'filled','SizeData', 80,'MarkerFaceColor', [0.5 0.5 0.5]);
alpha(0.7)

hold on
plot([-2.7 -1.3], [mean_cross(1) mean_cross(1)], 'Color', [0 0 0], 'MarkerSize',20, 'LineWidth', 3);
%plot([1.6 2.4], [mean_num(2) mean_num(2)] ,'Color', colors(1,:), 'MarkerSize',20, 'LineWidth', 3);

hold on
%   plot 30 textures with days

 
 category = {'_30_cross'};%, '30_iso'};

 for i=1:length(category)
     files=find_files( ['*' category{i} '*']);
     n = length(files);
     for j=1:n
         data(j) = load(files{j});
         nday(j)= size(data(j).correct_rate_total_sum,2);
         for k=1:size(data(j).correct_rate_total_sum,1)
             plot(data(j).correct_rate_total_sum(k,:),'-', 'Color', colors(i,:) ,'MarkerFaceColor', colors(i,:),'LineWidth', 2, 'MarkerSize',5);
         end% plot raw traces
         hold on
     end
     nday_min = min(nday);
     
     
     cr = [];
     for j=1:n
         cr = cat(1, cr, data(j).correct_rate_total_sum(:,1:nday_min));
     end
     cr_mean = mean(cr);
     cr_sem = std(cr)/sqrt(size(cr,1));
     %error_area(1:nday_min, cr_mean, cr_sem, min(colors(i,:)+0.25, [1 1 1]));
     hold on
     line(i) = plot(cr_mean,'-o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0],'LineWidth', 2);
     
     
     scatter( x*25, cross_7_test(:,2),'filled','SizeData', 80,'MarkerFaceColor', [0.5 0.5 0.5]);
     alpha(0.7)
     
     hold on
     plot([24.3 25.7], [mean_cross(2) mean_cross(2)], 'Color', [0 0 0], 'MarkerSize',20, 'LineWidth', 3);
     
     
     ylim([45 80]);
     xlim([-5 27+ 0.2]);
     
     % xlim([0.8 nday_min + 0.2]);
     set(gca,'XTick', 1:1:22 );
     %set(gca,'XTick', 1:1:14+0.2 );
     
     
     xlabel('Days');
     ylabel('Performance (% Correct)');
     % legend('mouse 1', 'mouse 2', 'mouse 3', 'mouse 4');
     set(gca, 'LineWidth', 2, 'FontSize', 14, 'TickDir','out');
     
     box off
     
 end
%%

% plot iso baseline training test
figure ( 'position',  [ 326         329        1130         420 ])
scatter( x*(-2), iso_7_test(:,1),'filled','SizeData', 80,'MarkerFaceColor', [0.5 0.5 0.5]);
alpha(0.7)

hold on
plot([-2.7 -1.3], [mean_iso(1) mean_iso(1)], 'Color', [0 0 0], 'MarkerSize',20, 'LineWidth', 3);
%plot([1.6 2.4], [mean_num(2) mean_num(2)] ,'Color', colors(1,:), 'MarkerSize',20, 'LineWidth', 3);

hold on
%   plot 30 textures with days

 
 category = {'_30_iso'};%, '30_iso'};

 for i=1:length(category)
     files=find_files( ['*' category{i} '*']);
     n = length(files);
     for j=1:n
         data(j) = load(files{j});
         nday(j)= size(data(j).correct_rate_total_sum,2);
         for k=1:size(data(j).correct_rate_total_sum,1)
             plot(data(j).correct_rate_total_sum(k,:),'-', 'Color', colors(i,:) ,'MarkerFaceColor', colors(i,:),'LineWidth', 2, 'MarkerSize',5);
         end% plot raw traces
         hold on
     end
     nday_min = min(nday);
     
     
     cr = [];
     for j=1:n
         cr = cat(1, cr, data(j).correct_rate_total_sum(:,1:nday_min));
     end
     cr_mean = mean(cr);
     cr_sem = std(cr)/sqrt(size(cr,1));
     %error_area(1:nday_min, cr_mean, cr_sem, min(colors(i,:)+0.25, [1 1 1]));
     hold on
     line(i) = plot(cr_mean,'-o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0],'LineWidth', 2);
     
     
     scatter( x*25, iso_7_test(:,2),'filled','SizeData', 80,'MarkerFaceColor', [0.5 0.5 0.5]);
     alpha(0.7)
     
     hold on
     plot([24.3 25.7], [mean_iso(2) mean_iso(2)], 'Color', [0 0 0], 'MarkerSize',20, 'LineWidth', 3);
     
     
     ylim([45 80]);
     xlim([-5 27+ 0.2]);
     
     % xlim([0.8 nday_min + 0.2]);
     set(gca,'XTick', 1:1:22 );
     %set(gca,'XTick', 1:1:14+0.2 );
     
     
     xlabel('Days');
     ylabel('Performance (% Correct)');
     % legend('mouse 1', 'mouse 2', 'mouse 3', 'mouse 4');
     set(gca, 'LineWidth', 2, 'FontSize', 14, 'TickDir','out');
     
     box off
     
 end







%% compare 30 cross vs. 30 iso
% 
% load('c_30_cross.mat')
% c_cross=ct_type_mean;

load('d_30_cross.mat')
d_cross=ct_type_mean;


load('e_30_cross.mat')
e_cross=ct_type_mean;

%c_cross;
cross_all = [ d_cross(1:2,:); e_cross];
%cross_all = [c_cross];
%cross_all = [ c_cross; d_cross(1:2,:); e_cross ];
%cross_all = [c_cross; d_cross; e_cross];%(1:2,:) d_cross(1:2,:)
mean_cross = mean(cross_all);
% 
% load('c_30_iso.mat')
% c_iso=ct_type_mean;

load('d_30_iso.mat')
e_iso=ct_type_mean;

load('e_30_iso.mat')
d_iso=ct_type_mean;

%c_iso;
iso_all = [ e_iso; d_iso(1:2,:)];
%iso_all = [c_iso;  d_iso(1:2,:); e_iso];
%iso_all = [c_iso];

%iso_all = [c_iso; e_iso; d_iso];
mean_iso = mean(iso_all);

figure ( 'position',  [  326   329   480   420]);
%plot(mean_iso, mean_cross,'o', 'Color', [1 1 1], 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',15,'LineWidth', 2);
scatter(mean_iso, mean_cross,'filled','SizeData', 150,'MarkerFaceColor', [0.5 0.5 0.5]);
alpha(0.5)

xlim([0.48 0.75]);
xlabel('Performances on 30 iso textures');
ylim([0.48 0.75]);
ylabel('Performances on 30 cross textures');

hold on
plot([0.48 0.73], [0.48 0.73], '--', 'Color', 'k', 'LineWidth', 2);
set(gca, 'LineWidth', 2, 'FontSize', 14,'TickDir','out');

%axis off
set(gca,'XTick', [0.5:0.1:0.7]);
set(gca,'YTick', [0.5:0.1:0.7]);
%axis equal
box off

%% summary of 30 texture task


% cross_7_test = xlsread('D:\experiment\ephys in mouse\paper_draft\behavior_Figure\30texture\7\cross_7_test.xlsx'); % before after cage c d e
% iso_7_test = xlsread('D:\experiment\ephys in mouse\paper_draft\behavior_Figure\30texture\7\iso_7_test.xlsx');
% %delete d 3,4    n=10
% cross_7_test=cross_7_test([1:4 7:end], :);
% iso_7_test=iso_7_test([1:4 7:end], :);


%delete cage c an cage d 3 4
cross_7_test = xlsread('D:\experiment\ephys in mouse\paper_draft\behavior_Figure\30texture\7\cross_7_test_delete_cage_c.xlsx'); % before after cage c d e
iso_7_test = xlsread('D:\experiment\ephys in mouse\paper_draft\behavior_Figure\30texture\7\iso_7_test_delete_cage_c.xlsx');


% c_30_cross_last_day=correct_rate_total_sum(:,end); 
% c_30_iso_last_day=correct_rate_total_sum(:,end); 
% 
% d_30_cross_last_day=correct_rate_total_sum(:,end); 
% d_30_iso_last_day=correct_rate_total_sum(:,end); 
% 
% e_30_cross_last_day=correct_rate_total_sum(:,end); 
% e_30_iso_last_day=correct_rate_total_sum(:,end); 

% cross_30_last_day= [c_30_cross_last_day; d_30_cross_last_day; e_30_cross_last_day];
% iso_30_last_day= [c_30_iso_last_day; d_30_iso_last_day;e_30_iso_last_day ];

load('D:\experiment\ephys in mouse\paper_draft\behavior_Figure\30texture\cross_iso_30_last_day.mat')

%delete cage c
cross_30_last_day = cross_30_last_day(5:end);
iso_30_last_day = iso_30_last_day(5:end);


%cross

figure ( 'position',  [  326   329   418   420]);  % 326   329   356   420
hold on
%plot(mean_iso, mean_cross,'o', 'Color', [1 1 1], 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',15,'LineWidth', 2);
plot(1, mean(cross_7_test(:,1)),'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',8,'LineWidth', 2);
errorbar(1, mean(cross_7_test(:,1)), std(cross_7_test(:,1))/sqrt(length(cross_7_test(:,1))), 'Color', [0.5 0.5 0.5],'LineWidth', 2 );

plot(1.5, mean(cross_30_last_day),'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',8,'LineWidth', 2);
errorbar(1.5, mean(cross_30_last_day), std(cross_30_last_day)/sqrt(length(cross_30_last_day)),'Color', [0.5 0.5 0.5],'LineWidth', 2 );

plot(2, mean(cross_7_test(:,2)),'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',8,'LineWidth', 2);
errorbar(2, mean(cross_7_test(:,2)), std(cross_7_test(:,2))/sqrt(length(cross_7_test(:,2))), 'Color', [0.5 0.5 0.5],'LineWidth', 2 );

ylim([45 80]);
set(gca, 'LineWidth', 2, 'FontSize', 14,'TickDir','out');


% iso
figure ( 'position',  [ 326   329   418   420]);
hold on
plot(1, mean(iso_7_test(:,1)),'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',8,'LineWidth', 2);
errorbar(1, mean(iso_7_test(:,1)), std(iso_7_test(:,1))/sqrt(length(iso_7_test(:,1))), 'Color', [0.5 0.5 0.5],'LineWidth', 2 );

plot(1.5, mean(iso_30_last_day),'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',8,'LineWidth', 2);
errorbar(1.5, mean(iso_30_last_day), std(iso_30_last_day)/sqrt(length(iso_30_last_day)),'Color', [0.5 0.5 0.5],'LineWidth', 2 );

plot(2, mean(iso_7_test(:,2)),'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',8,'LineWidth', 2);
errorbar(2, mean(iso_7_test(:,2)), std(iso_7_test(:,2))/sqrt(length(iso_7_test(:,2))), 'Color', [0.5 0.5 0.5],'LineWidth', 2 );

ylim([45 80]);
set(gca, 'LineWidth', 2, 'FontSize', 14,'TickDir','out');


box off

%%
% 
% cross_7_test = xlsread('D:\experiment\ephys in mouse\paper_draft\behavior_Figure\30texture\7\cross_7_test.xlsx'); % before after cage c d e
% iso_7_test = xlsread('D:\experiment\ephys in mouse\paper_draft\behavior_Figure\30texture\7\iso_7_test.xlsx');
% %delete d 3,4    n=10
cross_7_test = xlsread('D:\experiment\ephys in mouse\paper_draft\behavior_Figure\30texture\7\cross_7_test_delete_cage_c.xlsx'); % before after cage c d e
iso_7_test = xlsread('D:\experiment\ephys in mouse\paper_draft\behavior_Figure\30texture\7\iso_7_test_delete_cage_c.xlsx');


% cross_7_test=cross_7_test([1:4 7:end], :);
% iso_7_test=iso_7_test([1:4 7:end], :);


gray=[0.8 0.8 0.8];
symbor_size=8;
%x=ones(10,1);
 x=ones(6,1);

figure ( 'position',  [ 326   329    684   420 ])
%plot([1:2], all(:,1:2),'.-','Color',[0.5 0.5 0.5] ,'MarkerSize',20);%colors(1,:)
%plot( [1:2], all(:,1:2),'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2, 'MarkerSize',symbor_size);
scatter(x , cross_7_test(:,1),'filled','SizeData', 150,'MarkerFaceColor', [0.8 0.8 0.8]);
hold
scatter(x*2,  cross_7_test(:,2),'filled','SizeData', 150,'MarkerFaceColor', [0.8 0.8 0.8]);
plot([1:2], cross_7_test,'-','Color',[0.8 0.8 0.8] ,'MarkerSize',20,'LineWidth', 3);%colors(1,:)
%alpha(0.5)


hold on
plot([0.7 1.3], [mean(cross_7_test(:,1))  mean(cross_7_test(:,1)) ], 'Color', [0 0 0], 'LineWidth', 3);
plot([1.7 2.3], [mean(cross_7_test(:,2))  mean(cross_7_test(:,2)) ] ,'Color', [0 0 0], 'LineWidth', 3);

hold on 
%plot( [1:2], all(:,1:2),'-o', 'Color', gray, 'MarkerFaceColor', gray,'LineWidth', 2, 'MarkerSize',symbor_size);
scatter(x*4 , iso_7_test(:,1),'filled','SizeData', 150,'MarkerFaceColor', [0.8 0.8 0.8]);
hold on 
scatter(x*5,  iso_7_test(:,2),'filled','SizeData', 150,'MarkerFaceColor', [0.8 0.8 0.8]);
plot([4:5], iso_7_test,'-','Color',[0.8 0.8 0.8] ,'MarkerSize',20,'LineWidth', 3);%colors(1,:)
%alpha(0.5)


hold on
plot([3.7 4.3], [mean(iso_7_test(:,1))  mean(iso_7_test(:,1)) ], 'Color', [0 0 0], 'LineWidth', 3);
plot([4.7 5.3], [mean(iso_7_test(:,2))  mean(iso_7_test(:,2)) ] ,'Color', [0 0 0], 'LineWidth', 3);





xlim([0 6])
ylim([45 80]);
set(gca, 'LineWidth', 2, 'FontSize', 14,'TickDir','out');


box off

