#! /usr/bin/octave-cli -f 

% You should initialize with Init:Create_Input_Files.m , do the calculations with
% Calculate:Cap_Counts_AllPx10.m and analyze the results with Process:S2CAPout_AllPx10.m
% 'REFERENCE HERE'

% This program calculates the maximum number of GRBs in a given redshift range and sperical cap size 

% We think that this is optimized for speed - any major improvements are welcome! 

% INPUT DATA START

% Base parameters 
CAP_RANGE=0.02:0.01:1.80 ; % SPHERICAL CAP are divided by pi  
			   % !!WARNING: we impicitly assume everywhere that the step size is 0.01!! 

WINDOW_RANGE= 2:1:135; % WINDOW size in the redshift space

% Healpix grid centers N=6
% This file is created by Init:Create_Input_Files.m
healpix_file="HPIX/pixel_coords_map_nested_galactic_res6_xyz.bdat";

% Input file is created by Init:Create_Input_Files.m
% Input positional file is the command line arg(1) parameter. The file contains the (real or simulated) 
% grb_xyz array.  
% Serial number tells the type (0: original data, 1:1000 :simulations), not used here 
arg_list = argv(); grb_file=sscanf(arg_list{1},"%s"); 

% INPUT DATA END 

load(healpix_file); % read in of the healpix_xyz array
NHPIX=max(size(healpix_xyz)); % size of the Healpix grid 
load(grb_file); % read in of the grb_xyz array (ordered by redshift ! )

% Area selection integer: [galactic South, All, galactic North]
for area_selector=[-1 0 1]

	if(area_selector == 0 ) ; AREA="A"; grb=grb_xyz; 				 NGRB=max(size(grb)); end; %ALL data
	if(area_selector ==+1 ) ; AREA="N"; ix=find(grb_xyz(3,:)>=0); grb=grb_xyz(:,ix); NGRB=max(size(grb)); end; %Galactic NORTH
	if(area_selector ==-1 ) ; AREA="S"; ix=find(grb_xyz(3,:)<=0); grb=grb_xyz(:,ix); NGRB=max(size(grb)); end; %Galactic SOUTH

	% Here are the distances between the Healpix grid centers and the GRBs,
	% in degrees.  NGRB by NHPIX array 
	% Should be single precision for speed (real measurements' errors are
	% much bigger than these numerical errors) 
	angle =single(acos(grb'*healpix_xyz)*180/pi); 

	output_matrix=zeros(max(WINDOW_RANGE), max(floor(CAP_RANGE/0.01+0.5))); % Output array 
	% WINDOW selection: between 3-> NGRB-WINDOW+1 
	for WINDOW=WINDOW_RANGE
		% It is much faster to do the counting with the octave's highly
		% efficient matrix multiplication rather than do the cycling
		% and summing through the array's window sized ranges 

		% Counting matrix buildup:
		WINDOW_START_RANGE=1:NGRB-WINDOW+1 ; WN=max(size(WINDOW_START_RANGE));
		% This should be single precision for speed 
                window_selection_matrix=single(zeros(WN,NGRB));i=1;
		for WINDOW_START=WINDOW_START_RANGE
                        window_selection_matrix(i, WINDOW_START:WINDOW_START+WINDOW-1) = 1;
                        i++;
                end

		% Do the counting and maximum search
		for CAP=CAP_RANGE
			% CAP radius in degrees 
			CAP_RADIUS=acos( 1- CAP/2)* 180/pi;  
			% Which GRBs are inside the CAP_RADIUS at a given Healpix center? NGRB by NHPIX array of  0 or 1
			inside = single(angle<=CAP_RADIUS); 
			% Do the max search: we don't register the position (i.e. WINDOW_START)
			output_matrix(WINDOW,floor(CAP/0.01+0.5))=max( max(window_selection_matrix*inside)); 
		end
	end
	% save the results for a given area selector [galactic South, All, galactic North]
	filename=sprintf("S2CAPout/%s:k:%1s:AllPx10.dat", grb_file, AREA); save(filename,"-v7", "output_matrix");
end

