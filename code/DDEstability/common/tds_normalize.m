function tdsn=tds_normalize(tds)
%  tds_normalize: transform the delay-free E-matrix to identity
%
%     tdsn=tds_normalize(tds)
%
%  The function multiplies the state equation from the left such that
%  the matrix E corresponding to delay=0 is identity. If this is
%  already the case, no modifiction is done.
%
%
%  Example:
%    n=5;
%    A=randn(n); B=randn(n); C=randn(n); D=randn(n);
%    tds=tds_create({C,D},[0,1],{A,B},[1 3],'neutral');
%    tdsn=tds_normalize(tds);
%    shouldbezero=norm(tdsn.E{1}-eye(n))
%
tdsn=tds;
if (~isfield(tds,'E'))
    % If E is not specified, do nothing
    return;
end

I=find(tds.hE==0);
if (isempty(I))
    hE=tds.hE;
    error('No zero delay in hE found');
end

% Sum all the h=0 elements together
Esum=tds.E{I(1)};
for k=2:length(I)
    Esum=Esum+tds.E{I(k)};
end

% Do the multiplication
for k=1:length(tds.E)
    tdsn.E{k}=Esum\tdsn.E{k};
end
if (length(I)==1)
    % We know that this should be exactly eye
    if (issparse(tdsn.E{I}))
        tdsn.E{I}=speye(size(tdsn.E{I}));
    else
        tdsn.E{I}=eye(size(tdsn.E{I}));
    end
end

if (isfield(tdsn,'A'))
    for k=1:length(tds.A)
        tdsn.A{k}=Esum\tdsn.A{k};
    end
end

% Input matrices
if (isfield(tdsn,'B1'))
    for k=1:length(tdsn.B1)
        tdsn.B1{k}=Esum\tdsn.B1{k};
    end
end
if (isfield(tdsn,'B2'))
    for k=1:length(tdsn.B2)
        tdsn.B2{k}=Esum\tdsn.B2{k};
    end
end
