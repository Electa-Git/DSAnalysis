function tdso = tds_sort(tds)
% tds_sort: sort system matrices in ascending delay order.
% 
%    tds_compress(tds)

%  tdso is a struct with same fields as tds where the system delays are
%  sorted in ascending order and corresponding system matrices
%
%  Example:
%  tds=tds_create({3,2,1},[3 2 1]);
%  tdso = tds_sort(tds);

tdso=tds;
fields_mat={'E','A','B1','B2','C1','C2','D11','D12','D21','D22'};

for id=1:length(fields_mat)
    if isfield(tds,['h' fields_mat{id}])
        nhtmp=length(getfield(tds,['h' fields_mat{id}]));
        if (nhtmp>1)
            htmp=getfield(tds,['h' fields_mat{id}]);
            
            if ~issorted(htmp)
                [htmp_sorted,idx]=sort(htmp);
                mattmp=getfield(tds,fields_mat{id});
                tdso = setfield(tdso,fields_mat{id},mattmp(:,idx)); % rearrange matrices
                tdso = setfield(tdso,['h' fields_mat{id}],htmp_sorted); % sort delays                
            end
        end
    end
end