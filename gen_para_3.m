function [a] = gen_para_3(ang,num)
    temp=zeros(1,3);
    switch num
        case 1
            temp(3) = 0;
            temp(2) = 2*sin(ang)/sqrt(3);
            temp(1) = cos(ang) + sin(ang)/sqrt(3);
            
        case 3
            temp(1) = 0;
            temp(2) = sin(ang)/sqrt(3)-cos(ang);
            temp(3) = -2*cos(ang)-temp(2);
            
        case 2
            temp(2) = 0;
            temp(3) = -2*sin(ang)/sqrt(3);
            temp(1) = cos(ang) - temp(3)/sqrt(3);
        case 4
            temp(2) = 2/sqrt(3);
            temp(3) = -2*(sin(ang)-1)/sqrt(3);
            temp(1) = cos(ang) + temp(2)/2 +temp(3)/2;
        case 5
            temp(3) = 2/sqrt(3);
            temp(2) = 2*(sin(ang)+1)/sqrt(3);
            temp(1) = cos(ang) + temp(2)/2 +temp(3)/2;
        case 6
%             temp(1) = 1/2;
%             temp(2) = (2*sin(ang)/sqrt(3)-2*cos(ang)+1)/2;
%             temp(3) = temp(2)-2*sin(ang)/sqrt(3);
            temp(3) = 2/sqrt(3);
            temp(2) = 2*(sin(ang)+1)/sqrt(3);
            temp(1) = cos(ang) + temp(2)/2 +temp(3)/2;
            
        otherwise
            temp = 0;
    end

    temp_sum = sum(temp);
    a = temp/temp_sum;
end