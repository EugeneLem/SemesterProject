function X=create_X(x1,x2,x3,x4,u1,u2)
    %create a long list of all value where we want to identify the system

    l=length(x1)*length(x2)*length(x3)*length(x4)*length(u1)*length(u2);
    X=zeros(l,6);
    X(:,1)=kron(x1',ones(length(x2)*length(x3)*length(x4)*length(u1)*length(u2),1));
    X(:,2)=kron(ones(length(x1),1),kron(x2',ones(length(x3)*length(x4)*length(u1)*length(u2),1)));
    X(:,3)=kron(ones(length(x1)*length(x2),1),kron(x3',ones(length(x4)*length(u1)*length(u2),1)));
    X(:,4)=kron(ones(length(x1)*length(x2)*length(x3),1),kron(x4',ones(length(u1)*length(u2),1)));
    X(:,5)=kron(ones(length(x1)*length(x2)*length(x3)*length(x4),1),kron(u1',ones(length(u2),1)));
    X(:,6)=kron(ones(length(x1)*length(x2)*length(x3)*length(x4)*length(u1),1),u2');
end