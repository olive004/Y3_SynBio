%%
%step 1: analysis of abs600
%state which well of which sample should be discarded
%lb well 1 due to contamination
%c2 well 3 didnt grow well
%d2+ well 2 didnt grow well
clear all; close all;
ABS600=table2array(readtable('GROUP1.xlsx','Range','C7:Z88'));
ABS_LB=mean((ABS600(:,2:3)')); %treat this as zero, discard the first well data
ABS_A1=mean((ABS600(:,4:6)'))-ABS_LB; %value relative to LB, then take the average value 
ABS_B2=mean((ABS600(:,7:9)'))-ABS_LB; %B2_ABS=mean((B2_ABS'));
ABS_C2=mean((ABS600(:,10:11)'))-ABS_LB; %C2_ABS=mean((C2_ABS'));
ABS_D2=mean((ABS600(:,13:15)'))-ABS_LB; %D2_ABS=mean((D2_ABS'));
ABS_D2_PLUS=mean((ABS600(:,[16 18])'))-ABS_LB; %D2_PLUS_ABS=mean((D2_PLUS_ABS'));
ABS_S=mean((ABS600(:,19:21)'))-ABS_LB; %N_ABS=mean((N_ABS'));%standard expression strain 
ABS_N=mean((ABS600(:,22:24)'))-ABS_LB; %S_ABS=mean((S_ABS'));%no expression strain

figure(1); 
plot(0:15:1200,ABS_A1); 
hold all
plot(0:15:1200,ABS_B2); plot(0:15:1200,ABS_C2); plot(0:15:1200,ABS_D2); 
plot(0:15:1200,ABS_D2_PLUS); plot(0:15:1200,ABS_S); plot(0:15:1200,ABS_N);
xticks(0:60:1200);
xticklabels({'0','60','120','180','240','300','360','420','480','540','600','660','720','780','840','900','960','1020','1080',',1140','1200'});
xlabel('Time (mins)');
ylabel('ABS600 relative value (AU)');
legend('A1','B2','C2','D2','D2 PLUS','S','N');
title('ABS600 versus Time for A1, B2, C2, D2, D2+, S, and N promoters');

figure; plot(0:15:1200,ABS600(:,1)');
hold all
plot(0:15:1200,ABS600(:,2)');  plot(0:15:1200,ABS600(:,3)');
xticklabels({'0','60','120','180','240','300','360','420','480','540','600','660','720','780','840','900','960','1020','1080',',1140','1200'});
xlabel('Time (mins)');
xlabel('Time (mins)');
ylabel('ABS600 (AU)');
legend('well 1', 'well 2', 'well 3');
title('LB Control ABS600 Value')
%the exponential growth phase is between 120min and 180min
%for D2 plus and C2, use 480min to 540min
%%
%step 2: analysis of GFP
%state which well of which sample should be discarded
%lb well 1 due to contamination (too much GFP expression)
%c2 well 3 didnt grow properly and well 2 is contaminated
%d2+ well 2  didnt grow properly

gfp=table2array(readtable('GROUP1.xlsx','Range','B90:Z171'));
gfp_120_lb= mean(gfp(9,3:4)); %GFP of lb at 120min
gfp_180_lb= mean(gfp(13,3:4)); %GFP of lb at 180min
gfp_480_lb= mean(gfp(33,3:4));
gfp_540_lb= mean(gfp(37,3:4));
gfp_lb=gfp(:,2:4);
gfp_a1=mean((gfp(:,5:7)')); gfp_b2=mean((gfp(:,8:10)')); gfp_c2=mean((gfp(:,11:12)')); gfp_d2=mean((gfp(:,14:16)')); 
gfp_d2p=gfp(:,17)'; gfp_s=mean((gfp(:,20:22)')); gfp_n=mean((gfp(:,23:25)'));  
figure; plot(0:15:1200,gfp_lb(:,1)); 
hold on
plot(0:15:1200,gfp_lb(:,2)); plot(0:15:1200,gfp_lb(:,3))
xlabel('Time (mins)');
ylabel('GFP (AU)');
title('GFP for LB')
legend('well 1','well 2','well 3');

figure; plot(0:15:1200,gfp_a1);
hold on
plot(0:15:1200,gfp_b2); 
plot(0:15:1200,gfp_c2); 
plot(0:15:1200,gfp_d2); 
plot(0:15:1200,gfp_d2p); 
plot(0:15:1200,gfp_s); 
plot(0:15:1200,gfp_n);
xlabel('Time (mins)');
ylabel('GFP experssion (AU)');
title('GFP Versus Time For Different Samples');
legend({'a1','b2','c2','d2','d2p','s','n'});

%calculate GFP increment
gfp_120_a1= gfp(9,5:7)-gfp_120_lb; gfp_180_a1= gfp(13,5:7)-gfp_180_lb; gfp_inc_a1=gfp_180_a1-gfp_120_a1;%relative value, treat lb_gfp_120/180 as zero

gfp_120_b2= gfp(9,8:10)-gfp_120_lb; gfp_180_b2= gfp(13,8:10)-gfp_180_lb; gfp_inc_b2=gfp_180_b2-gfp_120_b2;

gfp_480_c2= gfp(33,13)-gfp_480_lb; gfp_540_c2= gfp(56,13)-gfp_540_lb; gfp_inc_c2=gfp_540_c2-gfp_480_c2;

gfp_120_d2= gfp(9,14:16)-gfp_120_lb; gfp_180_d2= gfp(16,14:16)-gfp_180_lb; gfp_inc_d2=gfp_180_d2-gfp_120_d2;

gfp_480_d2plus= gfp(33,19)-gfp_480_lb; gfp_540_d2plus= gfp(43,19)-gfp_540_lb; gfp_inc_d2plus=gfp_540_d2plus-gfp_480_d2plus; %d2+

gfp_120_s= gfp(9,20:22)-gfp_120_lb; gfp_180_s= gfp(13,20:22)-gfp_180_lb; gfp_inc_s=gfp_180_s-gfp_120_s;

gfp_120_n= gfp(9,23:25)-gfp_120_lb; gfp_180_n= gfp(13,23:25)-gfp_180_lb; gfp_inc_n=gfp_180_n-gfp_120_n;

%obtain abs300_t+30
abs_150_lb=mean(ABS600(11,2:3)); %abs600 of lb at 150min, ie ABS600_t+30
abs_510_lb=mean(ABS600(35,2:3));
abs_150_a1= mean(ABS600(11,4:6))-abs_150_lb; %abs600 of a1 at 150min
abs_150_b2= mean(ABS600(11,7:9))-abs_150_lb;
abs_510_c2= mean(ABS600(35,10:11))-abs_510_lb;
abs_150_d2= mean(ABS600(11,13:15))-abs_150_lb;
abs_510_d2plus= mean(ABS600(35,[16 18]))-abs_510_lb;
abs_150_s= mean(ABS600(11,19:21))-abs_150_lb;
abs_150_n= mean(ABS600(11,22:24))-abs_150_lb;

%Calibration information: ABS600*(0.5/0.58)=OD600
%calculate rate of GFP production per cell: GFP increment divided by
%ABS600_t+30
pro_rate_a1=gfp_inc_a1/(abs_150_a1*(0.5/0.58));
pro_rate_b2=gfp_inc_b2/(abs_150_b2*(0.5/0.58));
pro_rate_c2=gfp_inc_c2/(abs_510_c2*(0.5/0.58));
pro_rate_d2=gfp_inc_d2/(abs_150_d2*(0.5/0.58));
pro_rate_d2plus=gfp_inc_d2plus/(abs_510_d2plus*(0.5/0.58));
pro_rate_s=gfp_inc_s/(abs_150_s*(0.5/0.58));
pro_rate_n=gfp_inc_n/(abs_150_n*(0.5/0.58));


%calculate mean and standard error for d_gfp_p_cell, delta gfp per cell
%discard, the third well of c2 due to the E coli didnt grow properly
%the 2nd and 3rd of d2p due to the E coli didnt grow properly
mean_a1=mean(pro_rate_a1); st_error_a1=std(pro_rate_a1)/sqrt(3);
mean_b2=mean(pro_rate_b2); st_error_b2=std(pro_rate_b2)/sqrt(3);
mean_c2=mean(pro_rate_c2); st_error_c2=std(pro_rate_c2)/sqrt(1);
mean_d2=mean(pro_rate_d2); st_error_d2=std(pro_rate_d2)*sqrt(3);
mean_d2p=mean(pro_rate_d2plus(1)); st_error_d2p=pro_rate_d2plus(1);
mean_s=mean(pro_rate_s); st_error_s=std(pro_rate_s)/sqrt(3);
mean_n=mean(pro_rate_n); st_error_n=std(pro_rate_n)/sqrt(3);

%%
%step 3: determine the strength of the promoters in relative promoter unit
re_stren_a1=(mean_a1-mean_n)/(mean_s-mean_n); %relative strength of a1
re_stren_b2=(mean_b2-mean_n)/(mean_s-mean_n);
re_stren_c2=(mean_c2-mean_n)/(mean_s-mean_n);
re_stren_d2=(mean_d2-mean_n)/(mean_s-mean_n);
re_stren_d2p=(mean_d2p-mean_n)/(mean_s-mean_n);
re_stren_n=(mean_n-mean_n)/(mean_s-mean_n);
re_stren_s=(mean_s-mean_n)/(mean_s-mean_n);

name=categorical({'a1', 'b2', 'c2', 'd2', 'd2+' ,'n' ,'s'});
name=reordercats(name,{'a1', 'b2', 'c2', 'd2', 'd2+' ,'n' ,'s'});
value=[re_stren_a1 re_stren_b2 re_stren_c2 re_stren_d2 re_stren_d2p re_stren_n re_stren_s];
err=[0.1 0.05 0.01 0.001 0.03 0 0];
bar(name,value);
text(1:length(value),value,num2str(value'),'vert','bottom','horiz','center'); 
title('The Strength of the Promoters In Relative Promoter Unit');
hold on 
er = errorbar(name, value, err);
er.Color=[0 0 0];
er.LineStyle = 'none';
