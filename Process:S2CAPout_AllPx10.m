#! /usr/bin/octave-cli -f 

% You should initialize with Init:Create_Input_Files.m , do the calculations with
% Calculate:Cap_Counts_AllPx10.m and analyze the results with Process:S2CAPout_AllPx10.m
% You should reference the article: http://arxiv.org/abs/2312.10050

% This program aggrevates the Calculate:Cap_Counts_AllPx10.m output results from the S2CAPout/GRB+SIM/ directory 
% Data file should be like
% "S2CAPout/GRB+SIM/data-22-08-31_xyzpos:00000.bdat:k:A:AllPx10.dat",
% "S2CAPout/GRB+SIM/data-22-08-31_xyzpos:00001.bdat:k:A:AllPx10.dat", etc.
% Missing files are OK. 

% Base parameters: should be same as defined by Calculate:Cap_Counts_AllPx10.m  ! 

% SIMULATION_RANGE: it tells the type (0: original data, 1:1000 :simulations)
% could be larger than the # of simulations (i.e. missing simulations are OK) 
SIMULATION_RANGE=0:1000;

% Cap area sizes are divided by pi. 
% WARNING: we implicitly assume that the step size is 0.01 !! 
CAP_RANGE=0.02:0.01:1.20 ; 
NCAP=max(floor(CAP_RANGE./0.01+0.5) );

% Window size in the redshift space
WINDOW_RANGE= 2:1:120; 

% Area selection: [galactic South, All, galactic North]
for area_selector=[-1 0 1]

	% redshift szerinti radialis szanko + sapka 
	% input adatok:
	if(area_selector == 0 ) ; AREA="A"; end; %ALL data
        if(area_selector ==+1 ) ; AREA="N"; end; %Galactic NORTH
        if(area_selector ==-1 ) ; AREA="S"; end; %Galactic SOUTH

	% report array 
	report_array=zeros(max(size(SIMULATION_RANGE)), max(size(WINDOW_RANGE))*NCAP );

	for SIMNUM=SIMULATION_RANGE
		filename=sprintf("S2CAPout/GRB+SIM/data-22-08-31_xyzpos:%05d.bdat:k:%1s:AllPx10.dat",SIMNUM, AREA);
		% check the file: it could be missing because some cluster job error/power off, etc.
		if( exist (filename,"file") == 2) 
			load(filename);
			tmpdat=output_matrix(WINDOW_RANGE,:) ; 
			report_array(SIMNUM+1,:)=reshape((tmpdat),  1, max(size(WINDOW_RANGE))*NCAP );
		endif
	end

	% select only the valid outputs
	ix=find(sum(report_array')'); report_array=report_array(ix,:);

	% Calculate the MC number (cases with >= max value)
	out_number= reshape ( sum((sort(report_array(2:end,:),1)) >= repmat(report_array(1,:),size(report_array,1)-1,1) ,1),  max(size(WINDOW_RANGE)),NCAP  );
	% Calculate the MC rate  ((cases with >= max value)/(number of simulations))
	out_rate=out_number/(size(ix,1)-1);

	% Save it:
	filename=sprintf("S2CAPout/PPR2AllPx10:All::%1s:N:%05d.csv" , AREA, size(ix,1)); save(filename,'out_number');
	filename=sprintf("S2CAPout/PPR2AllPx10:Rate::%1s:N:%05d.csv", AREA, size(ix,1)); save(filename,  'out_rate');

	% Plot it:
	% ftitle=sprintf("cn:S2CAPout_AllPx10:All::%1s:N:%05d:", AREA,size(ix,1));
	% hf=figure (22); colormap('rainbow'); imagesc(log10(out_rate),'XData',[min(CAP_RANGE) max(CAP_RANGE)],'YData',[min(WINDOW_RANGE) max(WINDOW_RANGE)]); colorbar; title([ftitle "R2 log10 <= p0"]); xlabel('Cap Size'); ylabel('Window')
	% filename=sprintf("S2CAPout/PPR2AllPx10:Rate::%1s:N:%05dLog10.png",AREA,size(ix,1)); print (hf, filename , "-dpng");

end

