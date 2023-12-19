#! /usr/bin/octave-cli -f 

% You should initialize with Init:Create_Input_Files.m , do the calculations with
% Calculate:Cap_Counts_AllPx10.m and analyze the results with Process:S2CAPout_AllPx10.m
% You should reference the article: http://arxiv.org/abs/2312.10050

% This program is for preparing the input data 

% Healpix START
% A semiregular grid on the sphere is needed: Healpix is used here
% Calculate the Healpix point's xyz positions for a given resolution
for infix=[4 5 6 7] % Healpix resolution array - typically only one will be used 
	in_file =sprintf("HPIX/pixel_coords_map_nested_galactic_res%1d.dat",infix);
	out_file=sprintf("HPIX/pixel_coords_map_nested_galactic_res%1d_xyz.bdat",infix);
	healpix_array=load(in_file); % Healpix (index,l,b) read 
	healpix_l=healpix_array(:,2); healpix_b=healpix_array(:,3); healpix_index=healpix_array(:,1);
	healpix_x= cos(healpix_b/180*pi) .* cos(healpix_l /180*pi);
	healpix_y= cos(healpix_b/180*pi) .* sin(healpix_l /180*pi);
	healpix_z= sin(healpix_b/180*pi);
	% Merge the 3D vector 
	healpix_xyz=cat(2,healpix_x,healpix_y,healpix_z)'; 
	% Save the 3D coordinates of the Healpix centers
	save(out_file,'-binary', '-v7', "healpix_xyz");
end
% Healpix END


% GRB START
% Original GRB data
% We should sort the GRBs by redshift and write xyz positions according to that order
% real GRB data redshift+xyz positions on the unit sphere
grb_file='data-22-08-31::redshift:x:y:z.dat'
data=load(grb_file);
% sort by redshift 
[z ix]=sort(data(:,1)); 
% order by redshift
grb_x=data(ix,2);
grb_y=data(ix,2);
grb_z=data(ix,2);

% Joined GRB's xyz coordinates on the unit sphere, sorted by the redshift 
% I.e. radial information is in the order of the data! 
grb_xyz=cat(2,grb_x,grb_y,grb_z)'; 

% Real data is indexed with 0 and written out into the GRB+SIM directory
out_file=sprintf("GRB+SIM/%s_xyzpos:%05d.bdat", substr(grb_file,1, min (strfind(grb_file, ":"))-1),0 ); 
save(out_file,'-binary', '-v7', "grb_xyz");

% GRB END


% Mixing/MC Catalogues START
% For the MC calculactions we should do a permutation/mixing in the redshift
% (change the order of the data)

% The previous grb_xyz array (on the unit sphere) is the input
% Save the origial data 
original_grb_xyz=grb_xyz;
NGRB=max(size(original_grb_xyz));

% The mixing's index number will select between original (SIMUN=0) and/or the
% simulated data (SIMUNM>=0) 1000 simulations
for SIMNUM=1:1000
	[junk ix]=sort(rand(NGRB,1)); % we need a random order only - could be modified for a bootstrap 
	grb_xyz=original_grb_xyz(:,ix); % mixed xyz positions on the unit spheres - the radial 
					% information is in the order of the data! 
	% simulation is indexed with SIMNUM=1:1000, written into the GRB+SIM directory 
	out_file=sprintf("GRB+SIM/%s_xyzpos:%05d.bdat", substr(grb_file,1, min (strfind(grb_file, ":"))-1),SIMNUM ); 
	save(out_file,'-binary', '-v7', "grb_xyz");
end
% Mixing/MC Catalogues END

