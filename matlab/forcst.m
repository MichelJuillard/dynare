function [yf,var_yf]=forcst(dr,y0,k,var_list)
  global endo_nbr exo_nbr ykmin_ Sigma_e_ ex_ options_ lgy_
  
  options_.periods = k;
  make_ex_;
  yf = simult_(y0,dr,ex_(1:k,:),1);

  [A,B] = kalman_transition_matrix(dr_);
  
  sigma_u = B*Sigma_e_*B';
  sigma_y = 0;
  
  for i=1:k
    sigma_y = sigma_y+sigma_u;
    var_yf(i,dr.order_var) = diag(sigma_y(1:endo_nbr,1:endo_nbr))';
    if i == k
      break
    end
    sigma_u = A*sigma_u*A';
  end

  nvar = size(var_list,1);
  if nvar == 0
    nvar = endo_nbr;
    ivar = [1:nvar];
  else
    ivar=zeros(nvar,1);
    for i=1:nvar
      i_tmp = strmatch(var_list(i,:),lgy_,'exact');
      if isempty(i_tmp)
	disp(var_list(i,:));
	error (['One of the variable specified does not exist']) ;
      else
	ivar(i) = i_tmp;
      end
    end
  end

  for i=1:nvar
    my_subplot(i,nvar,2,3,'Forecasts');
    
    plot([-ykmin_+1:0],y0(ivar(i),1:ykmin_),'b-',...
	 [1:k],yf(ivar(i),ykmin_+1:end),'g-',...
	 [1:k],yf(ivar(i),ykmin_+1:end)'+2*sqrt(var_yf(:,ivar(i))),'g:',...
	 [1:k],yf(ivar(i),ykmin_+1:end)'-2*sqrt(var_yf(:,ivar(i))),'g:',...
	 [1 k],repmat(dr.ys(ivar(i)),1,2),'r-');
    title(lgy_(ivar(i),:));
  end










