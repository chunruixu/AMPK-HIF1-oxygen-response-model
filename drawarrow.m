function drawarrow(x,y,lineType,ax,string,color)
%drawarrow(xa,ya,'textarrow',gca,'Hypoxia','b');
%drawarrow(startpoint,endpoint,'String','color',linewidth);
% x, y: 箭头的起点和终点坐标，例如 [x_start, y_start; x_end, y_end]
% lineType: 线段类型，例如 'arrow' 表示画箭头
% ax: 坐标轴对象
% string: 在箭头旁边显示的文本
% color: 箭头和文本的颜色

% switch nargin
%     case 2
%         lineType = 'arrow';
%         ax = gca;
%     case 3
%         ax = gca;
% end

%调整坐标大小，ax = gca 返回当前图窗中的当前坐标区（或独立可视化）。使用 ax 获取和设置当前坐标区的属性
%当前坐标区x、 y 坐标轴范围
xlim = ax.XLim; ylim = ax.YLim;
xlimmin = xlim(1); xlimmax = xlim(2);
ylimmin = ylim(1); ylimmax = ylim(2);
%画的箭头x坐标在图传x坐标之外，把箭头x坐标设置为图窗最小x坐标
%xa = [Insert 35];%Insert=TimeOfRun=150
%调整箭头的起点终点坐标
if xlimmin>min(x(1),y(1)), xlimmin=min(x(1),y(1));end
if xlimmax<max(x(1),y(1)), xlimmax=max(x(1),y(1));end
if ylimmax>max(x(2),y(2)), ylimmax=max(x(2),y(2));end
if ylimmax<max(x(2),y(2)), ylimmax=max(x(2),y(2));end
ax.XLim = [xlimmin,xlimmax]; ax.YLim = [ylimmin,ylimmax];
xlim = ax.XLim; ylim = ax.YLim;
pos = ax.Position;
%pos(1) 是 Axes 对象左下角的 x 坐标。
%pos(2) 是 Axes 对象左下角的 y 坐标。
%pos(3) 是 Axes 对象的宽度。
%pos(4) 是 Axes 对象的高度
x_ratio = pos(3)/(xlim(2)-xlim(1));
%x_ratio 表示 Figure 中的一个单位长度（在 x 轴上）相当于 x 轴坐标轴上的多少长度。这个比例因子用于将实际的 x 轴坐标值映射到 MATLAB 中 Axes 对象的宽度上

y_ratio = pos(4)/(ylim(2)-ylim(1));%缩放比例
%-xlim(1)*x_ratio: 将数据坐标系中的 x 轴最小值 (xlim(1)) 转换为 Axes 对象坐标系中的 x 坐标
orig_pos = [-xlim(1)*x_ratio+pos(1),-ylim(1)*y_ratio+pos(2)];

%数值坐标*转换比例=图窗坐标
x=x.*[x_ratio,y_ratio]; y=y.*[x_ratio,y_ratio];
x=x+orig_pos; y=y+orig_pos;
annotation(lineType,[x(1),y(1)],[x(2),y(2)],'String',string,'FontSize',18,'Color',color)
end