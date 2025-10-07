%% Specific heat capacities of different SiO2 phases
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the thesis:
%
%Thanheiser, S.
%A Contribution to the Development of an Active Particle Thermal Energy 
%Storage System
%PhD Thesis, TU Wien, Austria, 2025
%
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.17288663
%
%
%
%This script creates Figure 4 in the thesis' framework showing specific
%heat capacities over temperature for the SiO2 phases quartz, tridymite,
%and cristobalite. 
%
%
%Required products, version 24.1:
%   - MATLAB
%   - Curve Fitting Toolbox
%Necessary files, classes, functions, and scripts:
%   - @SiO2


%% Set up figure
fig=figure(4);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


%% Plot lines
%Quartz
T=linspace(273.15,800,1e4);
plot(ax,T-273.15,SiO2.c_p(T,'quartz'),'Color',colors(1,:));

T=linspace(800,900,1e4);
plot(ax,T-273.15,SiO2.c_p(T,'quartz'),'Color',colors(1,:), ...
    'LineStyle',':','LineWidth',1);

T=linspace(900,2000,1e4);
plot(ax,T-273.15,SiO2.c_p(T,'quartz'),'Color',colors(1,:));


%Tridymite
T=linspace(273.15,370,1e4);
plot(ax,T-273.15,SiO2.c_p(T,'tridymite'),'Color',colors(2,:));

T=linspace(370,500,1e4);
plot(ax,T-273.15,SiO2.c_p(T,'tridymite'),'Color',colors(2,:), ...
    'LineStyle',':','LineWidth',1);

T=linspace(500,2000,1e4);
plot(ax,T-273.15,SiO2.c_p(T,'tridymite'),'Color',colors(2,:));


%Cristobalite
T=linspace(273.15,SiO2.TtransCri-50,1e4);
plot(ax,T-273.15,SiO2.c_p(T,'cristobalite'),'Color',colors(3,:));

T=linspace(SiO2.TtransCri-50,SiO2.TtransCri+50,1e4);
plot(ax,T-273.15,SiO2.c_p(T,'cristobalite'),'Color',colors(3,:), ...
    'LineStyle',':','LineWidth',1);

T=linspace(SiO2.TtransCri+50,2000,1e4);
plot(ax,T-273.15,SiO2.c_p(T,'cristobalite'),'Color',colors(3,:));


%Transition temperatures
T=SiO2.TtransQuartz;
plot(ax,repmat(T-273.15,1,2),[600,SiO2.c_p(T+1e-6,'quartz')],'Color','k');

T=SiO2.TtransTri;
plot(ax,repmat(T-273.15,1,2),[600,SiO2.c_p(T-1e-6,'tridymite')],'Color','k');

T=SiO2.TtransCri;
plot(ax,repmat(T-273.15,1,2),[600,SiO2.c_p(T+1e-6,'cristobalite')],'Color','k');


%% Set up legend
legItems=repmat(line(ax,'Visible','off'),4,1);
legItems(1)=plot(ax,NaN,NaN,'Color',colors(1,:));
legItems(2)=plot(ax,NaN,NaN,'Color',colors(2,:));
legItems(3)=plot(ax,NaN,NaN,'Color',colors(3,:));
legItems(4)=plot(ax,NaN,NaN,'Color','k','LineStyle',':','LineWidth',1);

lgd=legend(ax,legItems,...
    {'Quartz','Tridymite','Cristobalite','Extrapolation'},...
    'Location','best');

hold off


%% Text annotations for transition temperatures
text(ax,SiO2.TtransQuartz-273.15,1000,...
    {'\alpha \rightarrow \beta';'Quartz'},...
    'BackgroundColor','w',...
    'HorizontalAlignment','center');


text(ax,SiO2.TtransTri-273.15,680,...
    {'\alpha \rightarrow \beta';'Tridymite'},...
    'BackgroundColor','w',...
    'HorizontalAlignment','center');


text(ax,SiO2.TtransCri-273.15,850,...
    {'\alpha \rightarrow \beta';'Cristobalite'},...
    'BackgroundColor','w',...
    'HorizontalAlignment','center');


%% Axes configurations
%Parallel Kelvin axis
axTop=axes(t);
axTop.Color='none';
axTop.XAxisLocation='top';
axTop.YAxis.Visible='off';


%Axes limits
axTop.XLim=[273.15,1800];
ax.XLim=[273.15,1800]-273.15;


%Correct ticks
ax.XTick=0:200:1600;


%Axes labels
xlabel(axTop,'Temperature (K)');
xlabel(ax,'Temperature (°C)');
ylabel(ax,'Specific Heat Capacity (J/kgK)');


%% Export figure for manuscript
t.Units='centimeters';
t.OuterPosition=[0,0,17,8.5];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,'cpOverT.tiff','Resolution',600);




