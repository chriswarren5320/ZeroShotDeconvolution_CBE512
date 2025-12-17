function [U, V] = randomGridFlow(gridN, amp)
    U = (rand(gridN,gridN,'single')-0.5)*2*amp;
    V = (rand(gridN,gridN,'single')-0.5)*2*amp;
end