clc
close all
clear all

folder = '/home/dominic/Software/numericalAstrodynamicsTudatBundle/tudatBundle/tudat/Tudat/';

for i=1:6
          subplot(3,2,i)
   for j = 1:6
      data = load(strcat(folder,'closedLoopDopplerErrors_',num2str(i-1),'_',num2str(j-1),'.dat'));
      plot(data(:,1),data(:,2))
      hold on
      grid on
   end
end