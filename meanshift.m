function [ x ] = meanshift()

    mu = [0 0];
    S = [30,0;0 35];
    data = mvnrnd(mu,S,300);
    plot(data(:,1),data(:,2),'o');
    
    h = 2;
    up = [0,0];
    down = [0,0];
    x = [data(1,1),data(1,2)];
    
    hold on
    while 1
        plot(x(1),x(2),'r+');
        for i = 1:length(data)
            k = norm((x - data(i,:))/h).^2;
            up = up + data(i,:)*(1/sqrt(2*pi).*exp(-0.5*k));
            down = down +1/sqrt(2*pi).*exp(-0.5*k);
        end

        if(abs(x - up./down)<0.1)
            break;
        end
        x = up./down;
        up = 0;
        down = 0;
    end

end
