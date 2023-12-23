function [evps,VV1,H]=tds_arnoldi(tds,x0,N)
%  tds_arnoldi: Compute eigenvalues of delay system with cheb. Arnoldi
%
%  [evps,VV]=tds_arnoldi(tds,x0,N)
%  
%   Computes the eigenvalues of the time-delay system tds with the
%   delay Arnoldi scheme with starting vector x0 and N iterations.
%  
%   Example: 
% 
%   tds.A={rand(3),rand(3),rand(3)};
%   tds.hA=[0,1,2];
%   [evps,VV]=tds_arnoldi(tds,[1;1;1],20);
%   s=min(evps);
%   M=-s*eye(3)+tds.A{1}+tds.A{2}*exp(-s)+tds.A{3}*exp(-2*s);
%   shouldbezero=min(svd(M))
% 
% 
%
    
    
    n=length(tds.A{1});

    x0=x0/norm(x0); % normalize starting vector
    V=x0;
    H=[];

   
    taumax=max(tds.hA);
    
    %%  Matrix vector product with inverse of sum A_i
    Asum=tds.A{1}; 
    for k=2:length(tds.A)
        Asum=Asum+tds.A{k}; 
    end

    if (length(Asum)==1)
        Asuminvfun=@(x) x/Asum;
    else
        if (~issparse(Asum))
            Asuminvfun=@(x) Asum\x;
        else
             Asum=sparse(Asum);
            [L,U,P,Q]=lu(Asum);
            Asuminvfun=@(x) Q*(U\(L\(P*x)));
        end
    end
    
    %% Arnoldi iteration
    % fprintf([repmat('.', 1, ceil(N/5)) '\n'])
    fprintf('progress   this line''s duration    est. remaining time\n')
    startloop = tic;
    previous_duration = 0;
    total_lines = ceil(N/5);
    done_lines = 0;
    for k = 1:N
        if mod(k,5) == 0
            % fprintf('.')
            this_duration = toc(startloop);
            done_lines = done_lines + 1;
            % if previous_duration == 0
            %     previous_duration = this_duration;
            %     mean_increase_ratio = 1;
            % end
            % past_nb_lines = done_lines-1;
            remaining_nb_lines = total_lines-done_lines;
            % new_increase_ratio = this_duration/previous_duration;
            % if new_increase_ratio > 1.15
            %     new_increase_ratio = 1;
            % end
            % mean_increase_ratio = (mean_increase_ratio*past_nb_lines + new_increase_ratio)/(past_nb_lines+1);
            % remaining_duration = this_duration*mean_increase_ratio*(1-mean_increase_ratio^remaining_nb_lines)/(1-mean_increase_ratio);
            % if remaining_duration > 99999999
            %     remaining_duration = 99999999;
            % end
            remaining_duration = this_duration*remaining_nb_lines*1.1;
            % ---
            rd = remaining_duration;
            days  = floor(rd/(24*60*60));
            rd    = rd-days*24*60*60;
            hours = floor(rd/(60*60));
            rd    = rd-hours*60*60;
            mins  = floor(rd/60);
            rd    = rd-mins*60;
            if days > 0
                str_rd = [num2str(days) 'd ' num2str(hours) 'h ' num2str(mins) 'm ' num2str(floor(rd*100)/100) 's'];
            elseif hours > 0
                str_rd = [num2str(hours) 'h ' num2str(mins) 'm ' num2str(floor(rd*100)/100) 's'];
            elseif mins > 0
                str_rd = [num2str(mins) 'm ' num2str(floor(rd*100)/100) 's'];
            else
                str_rd = [num2str(floor(rd*100)/100) 's'];
            end
            % ---
            td = this_duration;
            days  = floor(td/(24*60*60));
            td    = td-days*24*60*60;
            hours = floor(td/(60*60));
            td    = td-hours*60*60;
            mins  = floor(td/60);
            td    = td-mins*60;
            if days > 0
                str_td = [num2str(days) 'd ' num2str(hours) 'h ' num2str(mins) 'm ' num2str(floor(td*100)/100) 's'];
            elseif hours > 0
                str_td = [num2str(hours) 'h ' num2str(mins) 'm ' num2str(floor(td*100)/100) 's'];
            elseif mins > 0
                str_td = [num2str(mins) 'm ' num2str(floor(td*100)/100) 's'];
            else
                str_td = [num2str(floor(td*100)/100) 's'];
            end
            fprintf(' %6.2f %%     %s       %s\n', k/N*100, str_td, str_rd)
            previous_duration = this_duration;
            startloop = tic;
        end    
        k0=size(V,1)/n; 
        
        % Apply matrix vector product to last column
        Y=reshape(V(:,end),n,k0);
        X=SigmaInvPiY(Y,tds,n,taumax,Asuminvfun);
        x=reshape(X,numel(X),1);
        
        % orthogonalization and H-update
        V=[V;zeros(n,size(V,2))];         
       
        [xorth,h]=gsr(V,x);
        H(1:k+1,k)=h;
        V=[V,xorth];              
    end
    % fprintf(' %6.2f %%  %6.2f s \n', k/N*100, toc(startloop))
    fprintf('\n')

    [vv,D]=eig(H(1:end-1,:));
    evps=1./diag(D);
    

    VV1=V(1:n,1:end-1)*vv;
    
   for k=1:size(VV1,2)
       VV1(:,k)=VV1(:,k)/norm(VV1(:,k));
   end

    
function [w,h]=gsr(v,w)
    h = v'*w ;
    w = w - v * h ;
    
    g = v'*w ;
    w = w - v * g ;
    h = h + g ;
    
    h = [h ; norm(w)] ;
    w = w / h(size(h,1)) ; 
end


function X=SigmaInvPiY(Yhat,tds,n,taumax,Asuminvfun)
%  The matrix vector product
%

    k=size(Yhat,2); 
    

    %  Construct the band matrix
    if (k==1)
        L=taumax*[0.5]; % Octave friendly 
    else
        dv=1:(k+1);
        L=diag(1./dv,0)+diag(-[1./dv(1:end-2)],2);
        L=L(1:end-1,:);
        L(1,1)=2;
        L=(taumax/4)*L;
        L=L(:,1:end-1);        
    end
    L=sparse(L);

    Z=Yhat*L';


    ysum=sum(Yhat,2);
    xx=ysum;
    for j=1:length(tds.hA)
        Tsum=zeros(size(xx));
        for i=1:size(Z,2)
            Tsum=Tsum+cheb_eval(1-2*tds.hA(j)/taumax,i)*Z(:,i);
        end
        xx=xx-tds.A{j}*Tsum;
    end
        
    xhat=Asuminvfun(xx);
    X=[xhat,Z];
end
    
function t=cheb_eval(x,n)    
%  Evaluate the chebyshev polynomial
    t=cos(n*acos(x));
end

end