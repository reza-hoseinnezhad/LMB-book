Model = gen_model;
%Truth = gen_truth(Model);
Meas=  gen_meas(Model,Truth);
Estimates =   run_filter(Model,Meas);
GeneratePlots(Truth,Meas,Estimates);