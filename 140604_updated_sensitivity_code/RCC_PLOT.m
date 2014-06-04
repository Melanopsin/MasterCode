%% Plot the ranked input (X - LHS matrix) vs the ranked output (Y) at column s (time point)
%% type: either linear-linear plot ('lin') or linear-log plot '(log')
%% RCC_var: labels of the parameters varied in the X (as legend)
%% y_var: label of the output variable
%% The Title of the plot is the Spearman correlation coefficient with the 
%% correspondent p-value and the time point chosen (s)
%% by Simeone Marino, June 5 2007 %%

function rcc_plot=RCC_PLOT(X,Y,s,type,RCC_var,y_var)
[A p]=corr(X,Y(s,:)','type','Spearman');
[a b]=size(X);
for i=1:b
    a=['[Spearman correlation , p-value] = ' '[' num2str(A(i)) ' , '  num2str(p(i)) '].'];% ' Time point=' num2str(s)]
    c=['ranking1(X(:,',num2str(i),'));'];
    %c2=['legend(var{i})'];
    if type=='log'
        figure,semilogy(eval(c),ranking1(Y(s,:)),'.'),title(a),legend(RCC_var{i}),...
            xlabel(RCC_var{i}),ylabel(y_var)
    

    elseif type=='lin'
        
        figure,plot(eval(c),ranking1(Y(s,:)),'ro', 'MarkerSize',5,'MarkerFaceColor', [1 0 0],'LineWidth',1.5);
        %red for I, blue for J
        %figure,plot(eval(c),ranking1(Y(s,:)),'bo', 'MarkerSize',5,'MarkerFaceColor', [0 0 1],'LineWidth',1.5);
                set(gca,'FontSize',20,'FontWeight','bold');
                set(gca, 'LineWidth', 3)
                
xlhand = get(gca,'xlabel');
set(xlhand,'string',RCC_var{i},'fontsize',20,'FontWeight','bold');
ylhand = get(gca,'ylabel');
set(ylhand,'string',y_var,'fontsize',20,'FontWeight','bold');
set(gca,'FontSize',20,'FontWeight','bold');
thand = get(gca,'title');
set(thand,'string',a,'fontsize',20,'FontWeight','bold');
set(gca,'FontSize',20,'FontWeight','bold');  

    end
end
