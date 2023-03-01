#include <cstdio>
#include <cmath>

#include "include/utils.h"
#include "include/enumerations.h"

int post_read_inputs()
{
	int err = 0;
	in_in_dir_name  = 'sample_inputs';
	in_out_dir_name = 'sample_outputs';

	in.R_max = 1.0e-4*100;
	in.tol_m = model_accu/10000;

	in.n_micro_elem = n_micro_elem;
	in.n_micro_step = n_micro_step;
	in.p_micro_step = p_micro_step;
	in.print_steps  = (in.p_micro_step-1:in.p_micro_step:in.n_micro_step-1);

	in.n_steps  = length(in.print_steps);
	lib.n_steps = in.n_steps;
	lib.stretch = (in.print_steps+1)*in.R_max/in.n_micro_step;

	in.raw_u1 = readmatrix([in_in_dir_name,'/sampling_results_1.csv']);
	in.raw_u2 = readmatrix([in_in_dir_name,'/sampling_results_2.csv']);
	in.raw_u3 = readmatrix([in_in_dir_name,'/sampling_results_3.csv']);

	in.raw_phi   = readmatrix([in_in_dir_name,'/sampling_results_2d_1.csv']);
	in.raw_theta = readmatrix([in_in_dir_name,'/sampling_results_2d_2.csv']);

	for k_cnt=1:length(in.print_steps)
	    k_ts = in.print_steps(k_cnt); %e.g.) 9 ~ 99

	    for i=1:n_sample
	        in.u1(n_sample*(k_cnt-1)+i,1) = in.raw_u1(i)/in.n_micro_step*(k_ts+1);
	        in.u2(n_sample*(k_cnt-1)+i,1) = in.raw_u2(i)/in.n_micro_step*(k_ts+1);
	        in.u3(n_sample*(k_cnt-1)+i,1) = in.raw_u3(i)/in.n_micro_step*(k_ts+1);

	        in.global_phi(n_sample*(k_cnt-1)+i,1)   = in.raw_phi(i);
	        in.global_theta(n_sample*(k_cnt-1)+i,1) = in.raw_theta(i);
	    end
	end
	fprintf('INFO : R_max = %e, n_sample=%d, n_elem=%d\n',in.R_max,n_sample,in.n_micro_elem);
	fprintf("INFO : accept when error is smaller than %2.2f\n",1-in.tol_m);

	return err;
}

int post_initialize_samples()
{
	int err = 0;

	max_out_error = zeros(n_sample*in.n_steps,1);
	out_error = zeros(n_sample*in.n_steps,1);
	out_error_micro = zeros(n_sample*in.n_steps,1);
	out_PP = zeros(n_sample*in.n_steps,1);
	out_PF = zeros(n_sample*in.n_steps,1);
	out_TP = zeros(n_sample*in.n_steps,1);
	out_TF = zeros(n_sample*in.n_steps,1);

	for i = 1:n_sample
	    fprintf('sample ID = %d\n',i);
	    for count_steps = 1:in.n_steps
	        p_fname = sprintf([in_out_dir_name,'/out_p',num2str(i),'/',num2str(np_micro),'CPU/pack_macro.out.',num2str(count_steps-1)]);
	        t_fname = sprintf([in_out_dir_name,'/out_t',num2str(i),'/',num2str(np_micro),'CPU/pack_macro.out.',num2str(count_steps-1)]);

	        p_fid=fopen(p_fname);
	        t_fid=fopen(t_fname);

	        p_temp_avg_data=fscanf(p_fid,'%f');
	        t_temp_avg_data=fscanf(t_fid,'%f');        

	        out_error_micro(n_sample*(count_steps-1)+i,1)=norm(t_temp_avg_data(25)-p_temp_avg_data(25))/norm(p_temp_avg_data(25));

	        // Microscale modeling error
	        // When small damage, no error (Use Taylor model)        
	        if t_temp_avg_data(25) <0.5 && p_temp_avg_data(25)<0.5
	            out_error_micro(n_sample*(count_steps-1)+i,1)=0.0;
	        end

	        // Macroscale modeling error
	        out_error(n_sample*(count_steps-1)+i,1)=norm(t_temp_avg_data(16:24)-p_temp_avg_data(16:24))./ norm(p_temp_avg_data(16:24));

	        // find the maximum error
	        max_out_error(n_sample*(count_steps-1)+i,1)=max(out_error(n_sample*(count_steps-1)+i,1),out_error_micro(n_sample*(count_steps-1)+i,1));

	        out_PP(n_sample*(count_steps-1)+i,1)=norm(p_temp_avg_data(16:24));
	        out_PF(n_sample*(count_steps-1)+i,1)=norm(p_temp_avg_data(1:9));
	        out_TP(n_sample*(count_steps-1)+i,1)=norm(t_temp_avg_data(16:24));
	        out_TF(n_sample*(count_steps-1)+i,1)=norm(t_temp_avg_data(1:9));

	        fclose(p_fid);
	        fclose(t_fid);
	    end // count_steps
	end // local_ID

	// data extraction from the largest pool with one-time conversion
	fprintf('INFO : sample_small_error directory generated\n')
	system('mkdir sample_small_error');
	fid=fopen(['./sample_small_error/error_',num2str(n_sample),'_pd.csv'],'w');
	fprintf(fid,'%1.12f\n',max_out_error);
	fclose(fid);
	clear out_error out_error_micro

	return err;

}

int post_determine_accept_reject_samples()
{
	int err = 0;

	count_id  = 1;
	ncount_id = 1;
	out_accept_id = -ones(1,1);
	out_reject_id = -ones(1,1);
	for count_steps = 1:in.n_steps
	    for i = 1:n_sample
	        if max_out_error(n_sample*(count_steps-1)+i) <= 1-in.tol_m            
	                out_accept_id(count_id) = n_sample*(count_steps-1)+i;
	                count_id = count_id+1;
	        else
	            out_reject_id(ncount_id) = n_sample*(count_steps-1)+i;
	            ncount_id = ncount_id+1;
	        end
	    end
	end
	clear temp_avg_PP temp_avg_TP;

	return err;

}


int post_generate_2d_layer()
{
	int err = 0;

	lib.point_r = zeros(lib.n_steps,1);
	lib.flag    = zeros(lib.n_steps,1);
	lib.n_a     = zeros(lib.n_steps,1);
	lib.B       = zeros(lib.n_steps,1);

	lib.AL = cell(lib.n_steps,1);
	lib.CL = cell(lib.n_steps,1);
	lib.XX = cell(lib.n_steps,1);
	lib.YY = cell(lib.n_steps,1);

	for stretch_count = 1:length(in.print_steps)
		map.a_count = 0;
	map.r_count = 0;
	map.a_id = zeros(1,1);
	map.r_id = zeros(1,1);

	for i = 1:n_sample
		in.local_phi(i)   = in.global_phi(n_sample*(stretch_count-1)+i);
	in.local_theta(i) = in.global_theta(n_sample*(stretch_count-1)+i);

	if max_out_error(n_sample*(stretch_count-1)+i) <= 1-in.tol_m
		map.a_count = map.a_count + 1;
	map.a_id(map.a_count) = i;
	else
		map.r_count = map.r_count + 1;
	map.r_id(map.r_count) = i;
	end
	end

	// PLOT: 2D-mapping
	gif.filename=['model_accuracy_map_SVM',num2str(in.tol_m*100),'_PmagDam.gif'];
	gif.h=figure(stretch_count+1);
	set(gif.h, 'Visible', 'off');

	// Decision boundary: SVM classifier
	if map.r_id(1)~=0 && length(map.r_id)>1 && map.a_id(1)~=0 && length(map.a_id)>1  // only if mixed cases
		if size(map.r_id)~=n_sample // strictly mixed cases (accepted-rejected)
			lib.flag(stretch_count) = 2; // mixed case 2
		map.X = [in.local_phi(map.r_id)' in.local_theta(map.r_id)';
			in.local_phi(map.a_id)' in.local_theta(map.a_id)'];

		map.class = [-ones(length(in.local_phi(map.r_id)'),1);
			ones(length(in.local_phi(map.a_id)'),1)];

		cl = fitcsvm(map.X,map.class,'KernelFunction','rbf','BoxConstraint',Inf,'ClassNames',[-1,1]);

		// save the data to global offline library
		lib.n_a(stretch_count)= length(cl.Alpha);
		lib.B(stretch_count)  = cl.Bias;
		lib.AL{stretch_count} = cl.Alpha;
		lib.CL{stretch_count} = map.class((cl.IsSupportVector));
		lib.XX{stretch_count} = map.X(cl.IsSupportVector,1);
		lib.YY{stretch_count} = map.X(cl.IsSupportVector,2);

		// Predict scores over the grid & Plot the data and the decision boundary
		d = 0.5*sqrt(2*pi*pi/n_sample); %h=1/2*(area/n_sample)^(1/2)

		[x1Grid,x2Grid] = meshgrid(min(map.X(:,1)):d:max(map.X(:,1)),...
			min(map.X(:,2)):d:max(map.X(:,2)));

		xGrid = [x1Grid(:),x2Grid(:)];

		[~,scores] = predict(cl,xGrid);
		map.class = [zeros(length(in.local_phi(map.r_id)'),1);
			ones(length(in.local_phi(map.a_id)'),1)];
		logict = logical(map.class);
		h(1) = plot(map.X(logict,1), map.X(logict,2),'r*');hold on; // accepted
		h(2) = plot(map.X(~logict,1),map.X(~logict,2),'b*');hold on; // rejected
		h(3) = plot(map.X(cl.IsSupportVector,1), map.X(cl.IsSupportVector,2),'ko'); hold on; // support vector
		contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),[0 0],'k'); hold off;            
		%legend(h,{'accepted','rejected','Support Vectors'}); %legend boxoff;
		legend(h,{'use Taylor','use PDE','Support Vectors'});
		end
		else
			if map.r_id(1)~=0 && length(map.r_id)>1 // Use PDE
				lib.flag(stretch_count) = 1; // Use PDE (rejected case 1)
			lib.n_a(stretch_count)  = 0;
			plot(in.local_phi(map.r_id),in.local_theta(map.r_id),'b*'); hold on;
			end
			if map.a_id(1)~=0 && length(map.a_id)>1 // Use Taylor
				lib.flag(stretch_count) = 0; // Use Taylor (accepted case 0)
			lib.n_a(stretch_count)  = 0;
			plot(in.local_phi(map.a_id),in.local_theta(map.a_id),'r*'); hold on;
			end
			end

			// additional plot setup
			xlim([0,pi*2]); ylim([0,pi]);
			set(gca,'FontSize',28,'FontName','Times');
			set(gcf,'color','w','position',[600, 600, 1600, 800]);
			pbaspect([2 1 1])
			make_gif(stretch_count,gif.filename,gif.h);
			close all;
			end

	return err;	
}


int post_write_offline_library()
{
	int err = 0;

	lib.fid = fopen('adaptive_params.in','w');

	fprintf(lib.fid,'# currently supported in develop_msnet branch\n');
	fprintf(lib.fid,'# read adaptive multiscale parameters\n\n');
	fprintf(lib.fid,'### -pre_micro option:\n');
	fprintf(lib.fid,'# thresholds setup: alpha_i * t_c and beta_j * u_c \n');
	fprintf(lib.fid,'# (alpha_0  alpha_1  alpha_2  alpha_3)\n');
	fprintf(lib.fid,'   0.0      0.0      0.0      0.0\n');
	fprintf(lib.fid,'# (beta_0  beta_1  beta_2  beta_3)\n');
	fprintf(lib.fid,'   0.0     0.0     0.0     0.0\n');
	fprintf(lib.fid,'# deactivatde when t_c or u_c is 0.0\n');
	fprintf(lib.fid,'# (t_c   u_c)\n');
	fprintf(lib.fid,'   0.0   0.0\n\n');
	fprintf(lib.fid,'### -post_micro option:\n');
	fprintf(lib.fid,'# required model accuracy \n');
	fprintf(lib.fid,'   %f',in.tol_m);
	fprintf(lib.fid,'\n');
	fprintf(lib.fid,'# the number of stretches: size of array\n');
	fprintf(lib.fid,'   %ld\n  ',lib.n_steps);
	for libI=1:lib.n_steps
		fprintf(lib.fid,' %2.8e',lib.stretch(libI));
	end
	fprintf(lib.fid,'\n');
	fprintf(lib.fid,'# the number of flags: size of array\n');
	fprintf(lib.fid,'# flag : 0 for Taylor, 1 for PDE, 2 requires parameters in [] \n');
	fprintf(lib.fid,'   %ld\n  ',lib.n_steps);
	lib.count_flag = 0;
	for libI=1:lib.n_steps
		fprintf(lib.fid,' %ld',lib.flag(libI));
	if lib.flag(libI)==2
		lib.count_flag = lib.count_flag + 1;
	end
	end
	fprintf(lib.fid,'\n');
	fprintf(lib.fid,'# the number of parameters: size of array\n');
	fprintf(lib.fid,'   %ld\n  ',lib.n_steps);
	for libI=1:lib.n_steps
		fprintf(lib.fid,' %ld',lib.n_a(libI));
	end
	fprintf(lib.fid,'\n');
	fprintf(lib.fid,'# Decision Boundary Parameters: flag-2 cases, size of array\n');
	fprintf(lib.fid,'# ( b  a  class  x1  x2 )\n');%flag 2-1-0 for mixed, rej, accept
	fprintf(lib.fid,'   %ld %ld\n',lib.count_flag,sum(lib.n_a*4)+lib.count_flag);
	if lib.count_flag > 0 // when MIXED cases
		for libI=1:lib.n_steps
			if lib.flag(libI) == 2
				fprintf(lib.fid,'   %2.8e ',lib.B(libI));
			fprintf(lib.fid,'%2.8e ',cell2mat(lib.AL(libI)));
			fprintf(lib.fid,'%2.8e ',cell2mat(lib.CL(libI)));
			fprintf(lib.fid,'%2.8e ',cell2mat(lib.XX(libI)));
			fprintf(lib.fid,'%2.8e ',cell2mat(lib.YY(libI)));
			fprintf(lib.fid,'\n');
			end
			end
			else // when NO mixed cases
			fprintf(lib.fid,' 0  0.0  0.0  0.0  0.0');
			fprintf(lib.fid,'\n');
			end
			fclose(lib.fid); // write library
			fprintf('\n...post_map_lib.m completed.\r\n');

	return err;
}