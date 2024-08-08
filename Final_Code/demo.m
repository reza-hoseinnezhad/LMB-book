Model = gen_model;
Truth = gen_truth(Model);
%plot_truth(Truth);
Meas=  gen_meas(Model,Truth);
[Estimates,LMBs] =   run_filter(Model,Meas);
%plot_truth_and_particles(Truth,LMBs);
GeneratePlots(Truth,Meas,Estimates);
