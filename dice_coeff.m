function similarity = dice_coeff(P, G, cluster_n)
%DICE_COEFF Computes Sørensen-Dice similarity coefficient for image
%           segmentation
%   similarity = dice_coeff(P, G) gets the proposed segmented image, P, and
%   the ground truth segmented image, G, as inputs and because we don't
%   know whether the labels are the same in both images or not, we should
%   consider all permutations of one of the images and compute dice
%   similarity coefficient for each. 

labels_p = 1 : cluster_n;
labels_g = unique(G);

p_perm = perms(labels_p);

max_dice = 0;
similarity = zeros(size(labels_p));
for l = 1 : size(p_perm, 1)
    G_new = zeros(size(G));
    for i = 1 : length(labels_g)
        G_new(G == labels_g(i)) = p_perm(l, i);
    end
    sim = dice(P, G_new);
    if sum(sqrt(sim(:).^2)) > max_dice
        max_dice = sum(sqrt(sim(:).^2));
        similarity = sim;
    end
end

end

