%% Plot X (input, LHS matrix) and Y (output) at column s (time point)
%% type: either linear-linear plot ('lin') or linear-log plot '(log')
%% var: labels of the parameters varied in the X (as legend)
%% The Title of the plot is the Pearson correlation coefficient
%% by Simeone Marino, June 5 2007 %%
function CC_PLOT(X,Y,s,type,PRCC_var,y_var)
[A p]=corr(X,Y(s,:)');
[a b]=size(X);
for i=1:b
    a=['[Pearson correlation , p-value] = ' '[' num2str(A(i)) ' , '  num2str(p(i)) '].'];% ' Time point=' num2str(s-1)]
    c=['X(:,',num2str(i),');'];
    



    if type=='lin'


figure,plot(eval(c),Y(s,:),'ro','MarkerSize',5,'MarkerFaceColor', [1 0 0],'LineWidth',1.5)% red for I
% figure,plot(eval(c),Y(s,:),'bo', 'MarkerSize',5,'MarkerFaceColor', [0 0 1],'LineWidth',1.5) %blue for J

set(gca, 'LineWidth', 3)


        set(gca,'FontSize',20,'FontWeight','bold');
xlhand = get(gca,'xlabel');
set(xlhand,'string',PRCC_var{i},'fontsize',20,'FontWeight','bold');
ylhand = get(gca,'ylabel');
set(ylhand,'string',y_var,'fontsize',20,'FontWeight','bold');
set(gca,'FontSize',20,'FontWeight','bold');
thand = get(gca,'title');
set(thand,'string',a,'fontsize',20,'FontWeight','bold');
set(gca,'FontSize',20,'FontWeight','bold');      
            
    elseif type=='log'
        figure,semilogy(eval(c),Y(s,:),'.'),title(a),legend(PRCC_var{i}),...%
            xlabel(PRCC_var{i}),ylabel(y_var);%eval(c6);
    end
end

