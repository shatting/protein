function info = plabel(rmse,sv,GDT)
if nargin == 3
    info = ['rmse= ',num2str(rmse),'  pot= ',num2str(sv),' GDT= ',num2str(GDT)]
else
    info = ['rmse= ',num2str(rmse),'  pot= ',num2str(sv)]
end