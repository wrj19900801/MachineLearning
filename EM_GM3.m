function X = EM_GM3( x,k )
% x now is a 2*1000 matrix
% k = 3 is the cluster we want to cluster
X=x'; %now X is a 1000*2 matrix
%X = [mvnrnd([1;1],eye(2),20);mvnrnd([10;10],eye(2),20);10*rand(20,2)];
[datasize dimension]=size(X);
mu = zeros(dimension,k)%mu is the center of GMM 2*3
sigma = zeros(dimension,dimension,k);
prior = zeros(1,k);%a 1*k matrix
%now we begin to intiliaze the parameters;
prior = repmat(1/k,1,k);
mu = max(max(X))*rand(dimension,k);
sigma = repmat(abs(max(max(X)))*eye(dimension),[1,1,k]);
mvnpdfs = zeros(datasize,k);

gamma = zeros(datasize,k);

for loop = 1:200
    mvnpdfs = zeros(datasize,k);        
    for i = 1:k
        mvnpdfs(:,i) = mvnpdf(X,mu(:,i)',sigma(:,:,i));
    end
    
    downsum = zeros(datasize,k);
    
    for ke = 1:k
        downsum(:,ke) = prior(ke)*mvnpdfs(:,ke);
    end
    
    for j = 1:datasize
        for ke=1:k
            gamma(j,ke) = prior(ke)*mvnpdfs(j,ke)/sum(downsum(j,:));
        end       
    end
 
    prior = sum(gamma)/datasize;
   
    for num = 1:k
        down = 0;
        up_sum=0;
        up_sum = repmat(gamma(:,num),1,2).*X;
        mu(:,num) = sum(up_sum)'/sum(gamma(:,num));
    end
    

    for num = 1:k
        up = zeros(2,2);
        down=0;
        for j=1:datasize
            up = up + gamma(j,num)*(X(j,:)-mu(:,num)')'*(X(j,:)-mu(:,num)');
        end
        sigma(:,:,num) = (up/sum(gamma(:,num)));
        sigma(:,:,num) = sigma(:,:,num)+.0001 * eye(dimension);
    end
end

hold on
for i = 1:datasize
    [max_temp,ind_temp] = max(gamma(i,:));
    if ind_temp == 1
        plot(X(i,1),X(i,2),'g.');
    end
    if ind_temp == 2
        plot(X(i,1),X(i,2),'b.');
    end
    if ind_temp == 3
        plot(X(i,1),X(i,2),'r.');
    end
end
hold off

end
