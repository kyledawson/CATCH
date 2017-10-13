% From Roger Stafford on Mathworks
% Finds closest value
%
% The following assumes that 'a' and 'b' are row vectors and that all the 
% elements of 'a' are finite. It obtains the two row vectors 'd' and 'ib' 
% of the same length as 'a'. Each element of 'd' is the absolute difference 
% between the corresponding element of 'a' and the nearest element of 'b'. 
% Each element of 'ib' is the index with respect to the 'b' vector of that 
% corresponding nearest 'b' element.
%
% modified by Kyle Dawson 02/14/2014
% added the 'n' argument, that sets the number of closest values.
% output is first closest, next closest, ..., n closest

function [ib,d] = closest(a,b,n)
 
a = reshape(a,1,length(a));
b = reshape(b,1,length(b));

if length(a)>length(b)
    b1 = a;
    b2 = b;
    a = b2;
    b = b1;
    clear b1 b2
end

if nargin<3
    n = 1;
end

ib = zeros([max(size(a)),n]);
d = ib;


for i1 = 1:n
 m = size(a,2); n = size(b,2);
 [~,p] = sort([a,b]);
 q = 1:m+n; q(p) = q;
 t = cumsum(p>m);
 r = 1:n; r(t(q(m+1:m+n))) = r;
 s = t(q(1:m));
 id = r(max(s,1));
 iu = r(min(s+1,n));
 [d(:,i1),it] = min([abs(a-b(id));abs(b(iu)-a)]);  % d is the abs distance
 ib(:,i1) = id+(it-1).*(iu-id);
 b(ib(:,i1))=inf;
 clearvars -except a b ib d
end