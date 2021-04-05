function [F] = myboxplot(boxes,names)
%Make boxplots of the given intervals
%   boxes is a 2 x n_box array that lists start and endpoints for each
%    interval
%   names is a 1 x n_box string cell that stores the names for the boxes
%   F is the handle for the figure generated

%% Initialize
whisker_length = .05;

n_box = size(boxes,2);
F = figure; hold on

xoffset = 1/n_box;

textshift = .05*(...
             max(max(boxes)) - min(min(boxes)) ...
            );
%% Plot the boxes
for i=1:n_box
    plot([i;i],boxes(:,i),'b','LineWidth',4,'Tag','Box');
    plot(i,mean(boxes(:,i)),'ko','MarkerSize',9,'MarkerFaceColor','w',...
                            'Tag','Middle_Circle');
    plot(i,mean(boxes(:,i)),'b.','MarkerSize',3,'MarkerFaceColor','b',...
                            'Tag','Middle_Dot');
    plot(i+whisker_length*[-1;1],repmat(boxes(1,i),2,1),'b','LineWidth',2,...
                            'Tag','Lower_Whisker');
    plot(i+whisker_length*[-1;1],repmat(boxes(2,i),2,1),'b','LineWidth',2,...
                            'Tag','Upper_Whisker');
end

for i=1:n_box
    textloc = F.CurrentAxes.YLim(1) - textshift;
    text(i,textloc,names{i},'HorizontalAlignment','center');
end

%% Finalize Figure Formatting
F.CurrentAxes.XLim(1) = xoffset;
F.CurrentAxes.XAxis.Visible = 'off';
