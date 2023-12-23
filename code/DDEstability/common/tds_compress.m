function tdso = tds_compress(tds)
% tds_compress: sums the system matrices in time-delay system with same 
% delays in the same type of system matrices.
%
%    tds_compress(tds)

%  tdso is a struct with same fields as tds where the system delays are
%  different in hE,hA,hB1,hC1,hD11,hB2,hC2,hD12,hD21,hD22 and the system
%  matrices with same delays in delay vectors are summed.
%
%  Example:
%  tds=tds_create({3,2,1},[3 2 2])
%  tdso = tds_compress(tds);

tdso=tds;
fields_mat={'E','A','B1','B2','C1','C2','D11','D12','D21','D22'};

for id=1:length(fields_mat)
    if isfield(tds,['h' fields_mat{id}])
        nhtmp=length(getfield(tds,['h' fields_mat{id}])); % number of delays in a field
        if (nhtmp>1) % more than one delay
            htmp=getfield(tds,['h' fields_mat{id}]); % get delay vector
            unihtmp=unique(htmp);
            nunihtmp=length(unihtmp);
            if (nhtmp~=nunihtmp) % there are same delays
                mattmp=getfield(tds,fields_mat{id});
                matouttmp={};
                tdso = rmfield(tdso,fields_mat{id}); % clear this system and delay fields in tdso
                tdso = rmfield(tdso,['h' fields_mat{id}]);
                for id1=1:nunihtmp
                   idx=find(htmp==unihtmp(id1)); 
                   if (length(idx)>1) % this delay appears more than one
                       smtmp=zeros(size(mattmp{1}));
                       for id2=1:length(idx) % add matrices with same delay
                           smtmp=smtmp+mattmp{idx(id2)};                       
                       end
                       matouttmp{id1}=smtmp; % set new system matrix
                   else
                       matouttmp{id1}=mattmp{idx}; % keep system matrix     
                   end
                end
                tdso = setfield(tdso,fields_mat{id},matouttmp); % set system field
                tdso = setfield(tdso,['h' fields_mat{id}],unihtmp); % set delay field
            end
        end
    end
end