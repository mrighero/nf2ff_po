% check dependencies
clc
clear all
close all

addpath('./auxiliary_routines_p');

% debug flag, now useless as it is hard coded to 0 (namely, no debug)
% within the routine 
debug_flag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set input parameters
fprintf(' --- START --- \n\n');
fprintf('Set parameters\n');
% mesh where sources are placed; sources are the v_n n=1, ... N mentioned
% before equation [7]
%data_in.source_mesh_name = '.\mesh\cylinder_arm_18GHz_l2.obj';
data_in.source_mesh_name = '.\mesh\box_arm_8GHz_l2.obj';
% mesh of the scatetring part; this is the tesselation T_scatterer
data_in.scatterer_mesh_name = '.\mesh\reflector_8GHz_8mm.obj';
% working frequency
data_in.frequency = 8e9;
% kind of measurement surface
data_in.measurement_kind = 'spherical';
data_in.planar = [];
% name of the file with the measured field samples; this is E^AUT(r_ell) in
% equation [2]
data_in.meas_field_name =  './field/SR40_SH4000_8d000GHz_CPQR_arm.mat';
% flag forcing the compuation of the far field (FF)
data_in.out_ff.compute = 1;
% minimum theta angle of the FF
data_in.out_ff.mintheta = 0;
% maximum theta angle of the FF
data_in.out_ff.maxtheta = 180;
% minimum phi angle of the FF
data_in.out_ff.minphi = 0;
% maximum theta angle of the FF
data_in.out_ff.maxphi = 360;
% thetaand phi angles step
data_in.out_ff.theta_step = 1;
data_in.out_ff.phi_step = 1;
% lsqr parameters: maximum number of iterations and target residue norm
data_in.lsqr.maxiter = 100;
data_in.lsqr.threshold = 1e-3;
% flag used during debug
data_in.spherical.force_step = 1;
data_in.process = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call to interpolate NF data
fprintf('nf2ff_po to compute FF...\n');
data_out = nf2ff_po(data_in, debug_flag);
fprintf('nf2ff_po to compute FF: DONE\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot of the FF, compared with the reference one

% reconstructed field
Field_rec = MR_Field_Format(data_out.out_ff.theta(:),...
    data_out.out_ff.phi(:),...
   [data_out.out_ff.Etheta(:), data_out.out_ff.Ephi(:)]);
% normalize to the maximum
Field_rec = MR_Field_Format(Field_rec.theta, Field_rec.phi,...
     [Field_rec.Etheta, Field_rec.Ephi]/max(Field_rec.Abs(:), [], 1));

% load reference FF
[theta_ref, phi_ref, Efield_ref] = MR_Read_FF2('./field/SR40_SH4000_8d000GHz.ffd');
Field_ref = MR_Field_Format(theta_ref, phi_ref, Efield_ref);
% normalize to the maximum
Field_ref = MR_Field_Format(Field_ref.theta, Field_ref.phi,...
     [Field_ref.Etheta, Field_ref.Ephi]/max(Field_ref.Abs(:), [], 1));

% plots
fprintf('Make FF plots:...\n');

% phi cuts where plots have to be made
CUTS = [0 45 90 135];

my_color = lines(2);

for kk = 1 : length(CUTS)
    
    h = figure('units','normalized','outerposition',[0.125 0 0.75 1]) ;
    hold all
    
    % reference
    index = find(abs(Field_ref.phi-CUTS(kk))<1e-1);
    a = Field_ref.theta(index);
    b = Field_ref.Y(index);
    c = Field_ref.X(index);

    index = find(abs(Field_ref.phi-(CUTS(kk)+180))<1e-1);
    a = [a(:); -Field_ref.theta(index(:))];
    b = [b(:); Field_ref.Y(index(:))];
    c = [c(:); Field_ref.X(index(:))];

    [a, index] = sort(a);

    plot(a, 20*log10(abs(b(index))), '-', 'color', my_color(1, :), 'linewidth', 2, 'displayname', 'CO ref')
    plot(a, 20*log10(abs(c(index))), '--', 'color', my_color(1, :), 'linewidth', 2, 'displayname', 'CX ref')
    
    % recconstructed
    index = find(abs(Field_rec.phi-CUTS(kk))<1e-1);
    a = Field_rec.theta(index);
    b = Field_rec.Y(index);
    c = Field_rec.X(index);
    
    index = find(abs(Field_rec.phi-(CUTS(kk)+180))<1e-1);
    a = [a(:); -Field_rec.theta(index(:))];
    b = [b(:); Field_rec.Y(index(:))];
    c = [c(:); Field_rec.X(index(:))];
    
    [a, index] = sort(a);
    
    plot(a, 20*log10(abs(b(index))), '-', 'color', my_color(2, :), 'linewidth', 2, 'displayname', 'CO rec')
    plot(a, 20*log10(abs(c(index))), '--', 'color', my_color(2, :), 'linewidth', 2, 'displayname', 'CX rec')
    
    % cosmetics
    xlabel('theta [deg]')
    ylabel('|E| [dBi]')
    title_string = sprintf('FF directivity pattern on the cut phi=%2.1f deg', CUTS(kk));
    title(title_string);
    grid on
    set(gca, 'fontsize', 20)
    legend('show')
       
    
end
fprintf('Make FF plots: DONE\n');
rmpath('./auxiliary_routines_p');


fprintf('\n --- END --- \n');
