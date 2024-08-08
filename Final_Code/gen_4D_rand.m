function Particles = gen_4D_rand(limits,npt)
    Particles = zeros(npt,4);
    xmin = limits(1);     xmax = limits(2);
    ymin = limits(3);     ymax = limits(4);
    xdotmin = limits(5);     xdotmax = limits(6);
    ydotmin = limits(7);     ydotmax = limits(8);
    Particles(:,1) = (xmax - xmin)*rand(npt,1) + xmin;
    Particles(:,2) = (ymax - ymin)*rand(npt,1) + ymin;
    Particles(:,3) = (xdotmax - xdotmin)*rand(npt,1) + xdotmin;
    Particles(:,4) = (ydotmax - ydotmin)*rand(npt,1) + ydotmin;
end