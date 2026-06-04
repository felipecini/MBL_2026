function [Start_ON,Stop_ON] = Detect_TTL_input(sg,threshold_ON)
% OD 01/05/2025
% Detect TTL onset/offset indices from analog signal sg using threshold
lgON = sg(:) > threshold_ON;
d = diff([false; lgON; false]);
Start_ON = find(d == 1);
Stop_ON = find(d == -1) - 1;
% fig
figure;
plot(-sg);hold on
bar(Start_ON,repmat(20,length(Start_ON),1),0.02,'r');
bar(Stop_ON,repmat(20,length(Start_ON),1),0.02,'g');

end