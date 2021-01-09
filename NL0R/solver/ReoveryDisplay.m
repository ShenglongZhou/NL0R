function ReoveryDisplay(xo,x,pos,ind)

    figure('Renderer', 'painters', 'Position', pos)
    stem(find(xo),xo(xo~=0),'bo-','MarkerSize',6, 'LineWidth',1),hold on
    stem(find(x),x(x~=0),'r*:', 'MarkerSize',4, 'LineWidth',1),hold on
    grid on, ymin= -0.1; ymax=0.2;
    xx  = [xo; x];
    if nnz(xx<0)>0, ymin= min(xx(xx<0))-0.1; end
    if nnz(xx>0)>0, ymax= max(xx(xx>0))+0.1; end   
    axis([1 length(x)  ymin  ymax])
    if ind
       st1   = strcat('Accuracy= ',num2str(norm(x-xo),2));
       wrong = nnz(xo)-nnz(find(x~=0 & xo~=0)); 
       st2   = strcat(', Sparsity= ',num2str(nnz(xo)),...
            ', Number of mis-indices= ',num2str(wrong)); 
       title(strcat(st1,st2))
       set(0,'DefaultAxesTitleFontWeight','normal');
       legend('Ground-Truth', 'Recovered', 'Location', 'best')
    end
end

