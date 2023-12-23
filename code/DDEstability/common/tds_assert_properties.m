function tds_assert_properties(tds,varargin)
% tds_assert_properties: checks that time-delay system has given properties
%
%   tds_assert_properties(tds,prop1,..., propk)
%
%  Checks if the time-delay system posseses the property specified
%  in the keywords prop1,..., propk. An error is thrown if not
%  all of them are fulfilled.
% 
%  Example:
%  tds=tds_create({rand(3)},0,{rand(3),rand(3)},[0 1],'neutral');
%  tds_assert_properties(tds,'lti','retarded'); % ok
%  tds_assert_properties(tds,'lti','neutral');  % ok
% 
%  tds=tds_create({rand(3),rand(3)},[0 1],{rand(3)},[0.1],'neutral');
%  tds_assert_properties(tds,'lti','retarded'); % error
% 
%  tds=tds_create({rand(3),rand(3)},[0 0],{rand(3)},[0.1],'neutral');
%  tds_assert_properties(tds,'lti','retarded','delay-free E'); % ok
%     
    

% For through all key words and check them 
    for i=1:length(varargin)
        isok=check_one_property(tds,varargin{i});
        if (isok<=0)
            if (isok==-1)
                error(['Unknown tds-property keyword: ',varargin{i}]);
            else
                error(['The function requires tds-property: ',varargin{i}]);
            end
        end
    end   
    % We're ok! 
end

function isok=check_one_property(tds,prop)
% Checks tds for property specified by string prop. returns 1 if
% ok, returns 0 if not fulfilled, returns -1 if keyword not identified
    switch prop
      case 'lti',  isok=1; % Currently always 
      
      case 'neutral',
        isok=0;
        E0=get_E0(tds);
        if (rank(E0)==length(E0))
            isok=1;
        end
      
      case 'retarded'
        if (~check_one_property(tds,'neutral'))
            isok=0;
            return;
        end
        isok=check_one_property(tds,'delay-free E');
      
      case 'delay-free IO'
        isok=check_one_property(tds,'delay-free B1') && ...
             check_one_property(tds,'delay-free B2') && ...
             check_one_property(tds,'delay-free C1') && ...
             check_one_property(tds,'delay-free C2') && ...
             check_one_property(tds,'delay-free D11') && ...
             check_one_property(tds,'delay-free D12') && ...
             check_one_property(tds,'delay-free D21') && ...
             check_one_property(tds,'delay-free D22');
      case 'SISO'
        
        maxion=0; % Compute sizes of all IOs 
        for k=1:length(tds.B1)
            if (length(tds.B1{k})>maxion)
               maxion=length(tds.B1{k}); 
            end
        end
        for k=1:length(tds.C1)
            if (length(tds.C1{k})>maxion)
                maxion=length(tds.C1{k});
            end
        end
        isok= (maxion<=1);
        
      otherwise
        isok=-1; %Unknown keyword
    end
    if (isok==-1) % Still unknown keyword 
        
        % Check the delay-free keywords
        if (1==findstr(prop,'delay-free '))
            k=regexp(prop, '\S*$');
            matname=prop(k:end); % Name of field for matrix
            if (~isfield(tds,matname))
                isok=-1; return; % Unknown matrix
            else
                hM=getfield(tds,['h',matname]);
                if (length(find(hM==0))==length(hM))
                    isok=1;
                else
                    isok=0;
                end
            end
            
        end
    end
end
 
    
% Returns the sum of E-matrices for the delay =0    
function E0=get_E0(tds)
    if (~isfield(tds,'E'))
        E0=eye(size(tds.A{1}));
        return;
    end
    if (length(tds.E)==0)
        E0=eye(size(tds.A{1}));
    else
        I=find(tds.hE==0);
        E0=tds.E{1};
        for k=2:length(I)
            E0=E0+tds.E{k};
        end
    end
end

       