function tds=tds_create(varargin)
% tds_create: creates time-delay object (tds) from input data
%   
%  call #1:
%  tds = tds_create(A,hA,B1,hB1,C1,hC1,D11,hD11,B2,hB2,C2,hC2,D12,hD12,D21,hD21,D22,hD22);
%  call #2:
%  tds = tds_create(E,hE,A,hA,B1,hB1,C1,hC1,D11,hD11,B2,hB2,C2,hC2,D12,hD12,D21,hD21,D22,hD22,'neutral');
%
%  where E,A,B1,C1,D11,B2,C2,D12,D21,D22 are cell arrays of system matrices
%  and hE,hA,hB1,hC1,hD11,hB2,hC2,hD12,hD21,hD22 are vectors of system
%  delays of the following state-space model
%
%  \sum_{i=1}^mE E{i}x'(t-hE{i}) = \sum_{i=1}^mA A{i} x(t-hA(i)) + \sum_{i=1}^mB1 B1{i} w(t-hB1(i)) + \sum_{i=1}^mB2 B2{i} u(t-hB2{i})
%  z(t) = \sum_{i=1}^mC1 C1{i} x(t-hC1(i)) + \sum_{i=1}^mD11 D11{i} w(t-hD11(i)) + \sum_{i=0}^nD12 D12{i} u(t-hD12(i))
%  y(t) = \sum_{i=1}^mC2 C2{i} x(t-hC2(i)) + \sum_{i=1}^mD21 D21{i} w(t-hD21(i)) + \sum_{i=0}^nD22 D22{i} u(t-hD22(i))
%
%  E and hE matrices are set when 'neutral' identifier is used at the end
%  of input parameters, otherwise E is set to identity matrix same size as
%  A matrices and hE is zero.
%
%  Examples:
%  % Create time-delay system without input or output: x'(t)=\sum_{i=1}^3 A{i}x(t-hA(i))
%  tds = tds_create({rand(3),rand(3), rand(3)},[0 1 2])
%
%  % Create time-delay system with input (u): x'(t)=\sum_{i=1}^3 A{i}x(t-hA(i))+ \sum_{i=1}^2 B1{i} u(t-hB1(i))
%  tds = tds_create({rand(3),rand(3), rand(3)},[0 1 2],{rand(3,1),rand(3,1)},[0 1.3])
%
%  % Create time-delay system with input (u) and output (y): 
%  %   x'(t)= \sum_{i=1}^3 A{i}x(t-hA(i))+ \sum_{i=1}^2 B1{i} u(t-hB1(i))
%  %   y(t) = \sum_{i=1}^2 C1{i} x(t-hC1(i)) + \sum_{i=1}^2 D11{i} u(t-hD11(i))
%  tds = tds_create({rand(3),rand(3), rand(3)},[0 1 2],{rand(3,1),rand(3,1)},[0 1.3],{rand(2,3),rand(2,3)},[0 1.1]...
%  ,{rand(2,1),rand(2,1)},[1 2])
%
%  % Create time-delay system in the form: E{1}x'(t)=\sum_{i=1}^3 A{i}x(t-hA(i))
%  tds = tds_create({2*eye(3)},0,{rand(3),rand(3), rand(3)},[0 1 2],'neutral')
%
%  % Create time-delay system with input in the form: \sum_{i=1}^2 E{i}x'(t-hE(i))=\sum_{i=1}^3 A{i}x(t-hA(i))+ \sum_{i=1}^2 B1{i} u(t-hB1(i))
%  tds = tds_create({2*eye(3),rand(3)},[0 1],{rand(3),rand(3), rand(3)},[0 1 2],{rand(3,1),rand(3,1)},[0 1.3],'neutral')
%
%  See also for examples: examples_tds_create.m
 
tds=[];
fields={'E','hE','A','hA','B1','hB1','C1','hC1','D11','hD11','B2','hB2','C2','hC2','D12','hD12','D21','hD21','D22','hD22'};

inputparam=varargin; % input parameters
nparam=length(inputparam);
idx=3; % index for field names, default no E and hE fields
if nparam==0
    error('No input parameters');
elseif ischar(inputparam{end}) && strcmpi(inputparam{end},'neutral')
    nparam=nparam-1; % exclude the last identifier 'neutral'
    inputparam=inputparam(1:nparam); % exclude last identifier
    idx=1; % the system has E and hE fields
end

% check whether there is missing matrix and fill it
id=1;
while (id<=nparam)
   if (id==nparam) && iscell(inputparam{nparam})
       nparam=nparam+1;
       if isempty(inputparam{nparam-1})
           inputparam{nparam} = []; % last delay is not entered
       else
           inputparam{nparam} = zeros(1,length(inputparam{nparam-1}));
       end
   elseif iscell(inputparam{id}) && iscell(inputparam{id+1})
       if isempty(inputparam(nparam))
           inputparam = {inputparam{1:id},[],inputparam{id+1:end}};
       else
           inputparam = {inputparam{1:id},zeros(1,length(inputparam{id})),inputparam{id+1:end}};
       end       
       nparam=nparam+1; % other delays are not given
   end
   id=id+1;
end

id=1; % input param
if (idx==3)
    tds = setfield(tds,'E',{}); % set matrix
    tds = setfield(tds,'hE',[]); % set matrix
end
while (idx<=length(fields))
    if (id<=nparam)
        tmpmat= inputparam{id};
        if ~iscell(tmpmat)
            error('System matrices in %s are not cell arrays',fields{idx});
        elseif ~isempty(tmpmat)
            for id1=1:length(tmpmat)
                if ~isnumeric(tmpmat{id1})
                    error('System matrix in %s is not numeric',fields{idx});
                end
            end
        end
        tds = setfield(tds,fields{idx},inputparam{id}); % set matrix
        
        tmpdel= inputparam{id+1};
        if ~isempty(tmpdel)
            if ~isvector(tmpdel)
                error('System delays in %s is not vector',fields{idx+1});
            elseif ~isnumeric(tmpdel)
                error('System delays in %s are not numeric',fields{idx+1});
            end        
        end
        tds = setfield(tds,fields{idx+1},inputparam{id+1}); % set delay        
    else
        tds = setfield(tds,fields{idx},{}); % set empty values
        tds = setfield(tds,fields{idx+1},[]);
    end
    
    idx=idx+2;
    id=id+2;
end

if isempty(tds.E)
    if isempty(tds.A)
        tds = setfield(tds,'E',{});
        tds = setfield(tds,'hE',[]);
    else
        tds = setfield(tds,'E',{eye(size(tds.A{1},2))});
        tds = setfield(tds,'hE',0);
    end
end