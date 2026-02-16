clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choices: 1 = Generate new; 2 = Load saved data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Choice = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate new ellipsoid %
%%%%%%%%%%%%%%%%%%%%%%%%%%

if Choice == 1
    
    X = [];
    Y = [];
    Z = [];

    eul = [pi/3 pi/6 pi/4]; % Euler Angles
    rotmZYX = eul2rotm(eul)

    xc = 25; % centre in mm
    yc = 25; % centre in mm
    zc = 25; % centre in mm

    a = 20; % radius in mm
    b = 15; % radius in mm
    c = 10; % radius in mm

    for i = 1:2:50
        for j = 1:2:50
            for k = 1:2:50
                if ((i-xc)^2/a^2 + (j-yc)^2/b^2 + (k-zc)^2/c^2 + randn/5) < 1
                    RotTrans = rotmZYX*[i-xc;j-yc;k-zc]+[xc;yc;zc]; 
                    X = [X; RotTrans(1)];
                    Y = [Y; RotTrans(2)];
                    Z = [Z; RotTrans(3)];
                end
            end
        end
    end

else
    load X.mat
    load Y.mat
    load Z.mat
end
    
%%%%%%%%%%%%%%%
% Plot Tumour %
%%%%%%%%%%%%%%%

% https://uk.mathworks.com/help/matlab/ref/convhull.html
[k2,av2] = convhull(X,Y,Z,'Simplify',true);
trisurf(k2,X,Y,Z,'FaceColor','cyan')
axis equal
xlim([0 50])
ylim([0 50])
zlim([0 50])
xlabel('x')
ylabel('y')
zlabel('z')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vertices for Students' Algorithm %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vertices = [X(k2(:,1)),Y(k2(:,1)),Z(k2(:,1))];
Vertices = [Vertices; X(k2(:,2)),Y(k2(:,2)),Z(k2(:,2))];
Vertices = [Vertices; X(k2(:,3)),Y(k2(:,3)),Z(k2(:,3))];
VerticesUnique = unique(Vertices,'rows');
hold on
scatter3(VerticesUnique(:,1),VerticesUnique(:,2),VerticesUnique(:,3),'r')
