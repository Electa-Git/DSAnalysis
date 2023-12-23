function tds_check_valid(tds)
% tds_metadata: checks that time-delay system is valid 
%
%   tds_metadata(tds)

%  Checks if time-delay system created by function tds_create satisfies the
%  following conditions:
%  - all fields for system matrices and delays are given
%  - system matrices are in a cell array with numerical content
%  - system delays are in a vector with numerical content
%  - system delays are nonnegative
%  - sizes of matrices inside a cell array is same
%  - the number of system delays are same as system matrices
%  - dimensions of E,A,B1,C1,D11,B2,C2,D12,D21,D22 are compatible
%
%  Example:
%  tds=tds_create({rand(2),rand(2)},[0 1]);
%  tds_check_valid(tds); % ok
%
%  tds=tds_create({rand(2),rand(2)},[0 -1]);
%  tds_check_valid(tds); % error
%
%  tds=tds_create({rand(2),rand(2)},[0 1 2]);
%  tds_check_valid(tds); % error

if nargin==0
    error('No input parameters');
elseif ~isstruct(tds)
    error('Time-delay system is not a struct');
end


%% check tds is LTI
sysType='lti'; % default

fields_mat={'E','A','B1','B2','C1','C2','D11','D12','D21','D22'};
fields_del=strcat('h',fields_mat);

% the system is LTI if there exists fields other than fields_mat
if isempty(fieldnames(tds))
    error('The struct has no fields');
else
    tdstmp=rmfield(tds,fields_mat(isfield(tds,fields_mat)));
    tdstmp=rmfield(tdstmp,fields_del(isfield(tds,fields_del)));
end

% if system is not LTI, throw error
if ~isempty(fieldnames(tdstmp))
    error('This program can handle LTI time-delay systems')
end % in future, update with other checks.

%% check LTI system is valid
if strcmp(sysType,'lti')
    % check all fields are defined
    for id=1:length(fields_mat)
        if ~isfield(tds,fields_mat{id})
            error('System matrix field %s of time-delay system is not defined',fields_mat{id});
        end
        if ~isfield(tds,fields_del{id})
            error('System delay field %s of time-delay system is not defined',fields_del{id});
        end
    end
    
    % fields of matrices and delays in tds
    fields_tdsmat=fields_mat(isfield(tds,fields_mat));
    fields_tdsdel=fields_del(isfield(tds,fields_del));
    
    % check each system matrix is cell array and its content is numeric
    % each system delay is vector and its content is numeric
    for id=1:length(fields_mat)
        tdsmat=getfield(tds,fields_mat{id});
        if ~iscell(tdsmat)
            error('The system matrices in field %s have to be a cell array',fields_mat{id});
        else
            for id1=1:length(tdsmat)
                tdsmat1=tdsmat{id1};
                if ~isnumeric(tdsmat1)
                    error('Matrix #%d in field %s is not numeric',id1,fields_mat{id});
                end
            end
        end
        
        tdsdel=getfield(tds,fields_del{id});
        if ~isempty(tdsdel)
            if ~isvector(tdsdel)
                error('The system delays in field %s have to be a vector',fields_del{id});
            else
                for id1=1:length(tdsdel)
                    tdsdel1=tdsdel(id1);
                    if ~isnumeric(tdsdel1)
                        error('Delay #%d in field %s is not numeric',id1,fields_del{id});
                    end
                end
            end
        end
    end
    
    mdtmp = struct;
    
    % sizes of matrices inside a cell array should be same
    for id=1:length(fields_tdsmat)
        tdsmat = getfield(tds,fields_tdsmat{id});
        
        mtdsmat=length(tdsmat);
        ntdsmat=[];
        if (mtdsmat==0)
            ntdsmat=[0 0];
        else                        
            for id1=1:mtdsmat
                ntdsmat(id1,:) = size(tdsmat{id1});
            end
        end
        
        if ~isempty(tdsmat)
            % check the dimensions of similar matrices are same
            if ~(((length(unique(ntdsmat(:,1)))==1) && (ntdsmat(1,1)==unique(ntdsmat(:,1)))) && ((length(unique(ntdsmat(:,2)))==1) && (ntdsmat(1,2)==unique(ntdsmat(:,2)))))
                error('The system matrices in %s do not have same dimensions',fields_tdsmat{id});
            elseif (strcmp(fields_tdsmat{id},'E') || strcmp(fields_tdsmat{id},'A')) && (ntdsmat(1,1)~=ntdsmat(1,2)) % check the matrices in A & E are square
                error('The system matrices in %s are not square matrices',fields_tdsmat{id});
            end
        end
        
        % dimension of the similar system matrices
        mdtmp = setfield(mdtmp,fields_tdsmat{id},ntdsmat(1,:));
        
        % check the corresponding delay vector is nonnegative
        tdsdel = getfield(tds,['h' fields_tdsmat{id}]);
        if any(tdsdel<0)
            error('a system delay of h%s is negative',fields_tdsmat{id});
        end
        
        % check the number of system delays are same as number of system matrices
        if (mtdsmat~=length(tdsdel))
            error('The number of matrices in %s and delays in %s are different',fields_tdsmat{id},fields_tdsdel{id});
        end
        % set the number of delay
        mdtmp = setfield(mdtmp,['h' fields_tdsmat{id}],mtdsmat);
    end
    
    % set the row and column dimension matrices
    rdimmat=[mdtmp.A(1) mdtmp.B1(1) mdtmp.B2(1);
        mdtmp.C1(1) mdtmp.D11(1) mdtmp.D12(1);
        mdtmp.C2(1) mdtmp.D21(1) mdtmp.D22(1)];
    
    cdimmat=[mdtmp.A(2) mdtmp.B1(2) mdtmp.B2(2);
        mdtmp.C1(2) mdtmp.D11(2) mdtmp.D12(2);
        mdtmp.C2(2) mdtmp.D21(2) mdtmp.D22(2)];
    
    dimmat_names={'A','B1','B2';'C1','D11','D12';'C2','D21','D22'};
    
    % dimension check between different matrices (except E)
    for idr=1:3
        for idc=1:3
            if (rdimmat(idr,idc)~=0)
                tmpr = rdimmat(idr,:);
                idnotcomp = find((tmpr~=rdimmat(idr,idc)) & (tmpr~=0), 1 );
                if ~isempty(idnotcomp)
                    error('The number of rows for %s and %s are not same',dimmat_names{idr,idc},dimmat_names{idr,idnotcomp});
                end
            end
            
            if (cdimmat(idr,idc)~=0)
                tmpc = cdimmat(:,idc);
                idnotcomp = find((tmpc~=cdimmat(idr,idc)) & (tmpc~=0), 1 );
                if ~isempty(idnotcomp)
                    error('The number of columns for %s and %s are not same',dimmat_names{idr,idc},dimmat_names{idnotcomp,idc});
                end
            end
        end
    end
    
    % dimension check for E
    if (mdtmp.hE~=0)
        idnotcomp = find(mdtmp.E(1,1)~=rdimmat(1,:) & (rdimmat(1,:)~=0),1);
        if ~isempty(idnotcomp)
            error('The number of rows for E and %s are not same',dimmat_names{1,idnotcomp});
        end
        
        idnotcomp = find(mdtmp.E(1,2)~=cdimmat(:,1) & (cdimmat(:,1)~=0),1);
        if ~isempty(idnotcomp)
            error('The number of columns for E and %s are not same',dimmat_names{idnotcomp,1});
        end
    elseif ~isempty(tds.A)
        error('Set E matrix to identity with same size as A matrices and hE to 0');
    end
end
