function inverted = invertLU(A)
    
    s = size(A);
    s = s(1);

    % LU factorization
    [L, U, P] = lu(A);
    
%     % creating identity matrix
%     I = eye(s);
%     
%     inverted = zeros(s);
%     X = zeros(s);
% 
%     
%     for k=1:s
%         X(1,k) = I(1,k);
%         for m=2:s
%             X(m,k) = (I(m,k)-L(m,1:m-1)*X(1:m-1,k));
%         end
%         
%         inverted(s,k) = X(s,k)./U(s,s);
%         
%         for m=s-1:-1:1
%             inverted(m,k) = (X(m,k)-U(m,m+1:s)*inverted(m+1:s,k))./U(m,m);
%         end
%     end
%    

I=eye(size(A));
s=size(A,1);
inverted=zeros(size(A));
for i=1:s
    b=I(:,i);
    inverted(:,i)=TriangleBackwardSub(U,TriangleForwardSub(L,P*b));
end

    
end

