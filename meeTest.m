clc
clear all
close all

directory = '/home/dominic/Software/numericalAstrodynamicsTudatBundle/tudatBundle/tudat/Tudat/';

variablesDiff=cell(5,1);
for i=1:5
    
    orbit = load(strcat(directory,'meeThrustState_',num2str(i-1),'.dat'));
    variables = load(strcat(directory,'meeThrustDep_',num2str(i-1),'.dat'));
    variablesDiff{i} = zeros(size(variables));
    
    for j=1:max(size(variables))
        variablesDiff{i}(j,:)=variables(j,:)-variables(1,:);
    end
    
     figure(1)
     subplot(2,3,i)
     plot3(orbit(:,2),orbit(:,3),orbit(:,4))
     axis equal
     
     figure(2)
     for j=1:5
        subplot(3,2,j)
        semilogy(orbit(:,1),abs(variablesDiff{i}(:,j+1)))
        hold on
        grid on
     end
    
     figure(3)
     for j=1:5
        subplot(3,2,j)
        plot(orbit(:,1),(variablesDiff{i}(:,j+1)))
        hold on
        grid on
    end
%     
%     figure(i+6)    
%     plot(orbit(:,1),dot(orbit(:,5:7)',variables(:,8:10)')')
% 
%     figure(i+12)  
%     
%     subplot(1,2,1)    
%     for j=1:3
%         plot(orbit(:,1),orbit(:,j+1))
%         hold on
%     end
%     
%     subplot(1,2,2)    
%     for j=1:3
%         plot(orbit(:,1),variables(:,j+7))
%         hold on
%     end    
end