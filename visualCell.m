%visualize the supercell and force
clf
hold on

for inc=1:ndl
    for ia=1:4
        scatter3(RR{ia+(inc-1)*8}(1),RR{ia+(inc-1)*8}(2),RR{ia+(inc-1)*8}(3),'fill','r');
    end
    
    for ia=5:8
        scatter3(RR{ia+(inc-1)*8}(1),RR{ia+(inc-1)*8}(2),RR{ia+(inc-1)*8}(3),'fill','b')
    end
end

% for ia=1:na
%     quiver3(RR{ia}(1),RR{ia}(2),RR{ia}(3),Force{ia}(1),Force{ia}(2),Force{ia}(3),'g')
% end

daspect([1 1 1])