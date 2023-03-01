#include <cstdio>
#include <cmath>

#include "include/utils.h"
#include "include/enumerations.h"

// generate samples based on "jump displacements"

int pre_build(PreVars *prv)
{
	int err = 0;
	double *print_time_step;


	// user parameters, PGFem3D inputs
	char file_name = 'sample_';
	char file_base = 'pack';

	double R_max = 0.005; // 5 percent of 0.1mm
	double dt = 2.5e-4; // rate 2
	int n_micro_step = 200;
	int print_by_n = prv->p_micro_step;
	double dudt= (R_max/n_micro_step)/dt;

	// start
	switch (prv->opt_ts_print) {
	    case 0 // all steps
	    	print_time_step = calloc(double,n_micro_step);
	        print_time_step = 0:1:n_micro_step-1;
	        break;
	    case 1 // initial & last
	        print_time_step = [0 n_micro_step-1];
	        break;
	    case 2 // 10 steps + last
	        print_time_step = [print_by_n-1:print_by_n:n_micro_step-1]; 
	        break;
	    default:
	        print_time_step = 0:1:n_micro_step-1;
	        break;
	}

	//<<-- START INPUT COMPUTING
	if exist([file_name,'1'],'dir')
	    fprintf('...cleaning sample_* directories \n');
	    system('rm -r ./',file_name,'* ');
	    fprintf('...complete \n');
	end

	// STEP 1: samples from uniform distribution u~U(,)
	norm_result = zeros(1,n_sample);
	raw_result = zeros(3,n_sample);  // raw data

	raw_plot = zeros(size(raw_result,1), size(raw_result,2)); // time step
	raw_plot_2d= zeros(2,n_sample); // raw data of phi, theta

	// quasi-random parameter
	rng default  // For reproducibility
	p = haltonset(2,'Skip',1e3,'Leap',1e2);
	p = scramble(p,'RR2');
	X0 = net(p,n_sample); %X0 = p(1:500,:);


    err += post_generate_samples();
    
	// tensile sampling
	j=1;
	n_trials=1;
	while (j < n_sample +1) {		

	    switch opt_R_max // r sampling
	        case 0 // surface
	            R = R_max;
	            du = R/n_micro_step;
	            fprintf('...generate %d/%d samples, R_max =%e, %d steps, ||du||=%e, dt=%e, ||du||/dt=%e\r\n',...
	            j,n_sample,R,n_micro_step,du,dt,dudt);
	        case 1 // radius
	            R = unifrnd(0,R_max);
	            du = R/n_micro_step;
	            dt = du/dudt;
	            fprintf('...generate %d/%d samples, R_max =%e, %d steps, ||du||=%e, dt=%e, ||du||/dt=%e\r\n',...
	            j,n_sample,R,n_micro_step,du,dt,dudt);
	        case 2 // volume
	            R = R_max * (unifrnd(0,1))^(1/3);
	            du = R/n_micro_step;
	            dt = du/dudt;            
	            fprintf('...generate %d/%d samples, R_max =%e, %d steps, ||du||=%e, dt=%e, ||du||/dt=%e\r\n',...
	            j,n_sample,R,n_micro_step,du,dt,dudt);
	        otherwise
	            fprintf('ERROR: please specify the option for R\n');
	            return
	    end

	    
	    //phi = unifrnd(0,2*pi);
	    //theta = acos(1-unifrnd(0,2));
	    //phi = pi/4*X0(j,1);
	    phi = 2*pi*X0(j,1);
	    theta = pi*X0(j,2);
	    
	    u1=R*cos(phi)*sin(theta);
	    u2=R*sin(phi)*sin(theta);
	    u3=R*cos(theta);
	    
	    norm_xx = sqrt(u1^2+u2^2+u3^2);
	    
	    // if norm_xx <= R_max
	        raw_result(1,j)=u1;
	        raw_result(2,j)=u2;
	        raw_result(3,j)=u3;        
	        
	        // norm
	        norm_result(1,j)= norm_xx;        
	        raw_plot(1,j)=u1;
	        raw_plot(2,j)=u2;
	        raw_plot(3,j)=u3;

	        raw_plot_2d(1,j)=phi;
	        raw_plot_2d(2,j)=theta;
	    
	        j=j+1;
	    %end
	    n_trials = n_trials +1;
	end
	fprintf("...complete. %d samples from %d trials \r\n",n_sample,n_trials);


	fidd = fopen('sampling_results_1.csv','w');
	fprintf(fidd,'%f\n',raw_plot(1,:));
	fclose(fidd);
	fidd = fopen('sampling_results_2.csv','w');
	fprintf(fidd,'%f\n',raw_plot(2,:));
	fclose(fidd);
	fidd = fopen('sampling_results_3.csv','w');
	fprintf(fidd,'%f\n',raw_plot(3,:));
	fclose(fidd);

	fidd = fopen('sampling_results_2d_1.csv','w');
	fprintf(fidd,'%f\n',raw_plot_2d(1,:));
	fclose(fidd);
	fidd = fopen('sampling_results_2d_2.csv','w');
	fprintf(fidd,'%f\n',raw_plot_2d(2,:));
	fclose(fidd);



	// STEP 2 :: generate PGFem3D input file
	err += pre_generate_pgfem_input();

	
	fprintf('...note: make sure numerical_param.json for details\n');
	fprintf('...complete. PGFem3D input type (%d) generated\n',opt_par);

	} // end of while j< ... loop

	return err;

}


int pre_generate_pgfem_input();

{
	int err =0;

	for (int sample_id = 1; sample_id < n_sample+1; sample_id++) {
		    
		    dir_name = sprintf('sample_%d',sample_id);
		    if (~exist(dir_name, 'dir')){
		        mkdir(dir_name)
		    }
		    
		    if (opt_par == 1){
		        if (~exist([dir_name,'/in_st_p'], 'dir')){
		            mkdir([dir_name,'/in_st_p']);
		        }
		        if (~exist([dir_name,'/in_st_t'], 'dir')){
		            mkdir([dir_name,'/in_st_t']);
		        }
		    }
		    
		    d_jump = 1/n_micro_step*[raw_result(1,sample_id);
		                             raw_result(2,sample_id);
		                             raw_result(3,sample_id)];
		    
		    cnt = 1:n_micro_step;
		    total_u= zeros(3,n_micro_step);
		    total_t= zeros(1,n_micro_step);
		    
		    for i = 1:n_micro_step
		        total_u(1,i) = d_jump(1)*cnt(i);
		        total_u(2,i) = d_jump(2)*cnt(i);
		        total_u(3,i) = d_jump(3)*cnt(i);
		        total_t(i)   = dt*cnt(i);
		    end
		    
		    switch opt_par // inside sample_id loop
		        case 0 // human-readable inputs
		            // Mechanical loads, load file name
		            fid = fopen([dir_name,'/Mechanical_load.json'], 'w');
		            
		            // write loadsteps, sizes
		            fprintf(fid, '{\n');
		            fprintf(fid, '  "number_of_loading_steps": %d, \n', n_micro_step-1);
		            fprintf(fid, '  "#": "the next 2 are ignored if number_of_loading_steps is 0",\n');
		            fprintf(fid, '  "loading_time_steps": \n');
		            fprintf(fid, '   [');
		            fprintf(fid, '%d, ',1:1:n_micro_step-2);
		            fprintf(fid, '%d ',n_micro_step-1);
		            fprintf(fid, '], \n');
		            
		            // start increment
		            fprintf(fid, '  "loading_increments":\n  [');
		            fprintf(fid,'\n   ['); fprintf(fid, '%e, ',total_u(1,2:end-1)); fprintf(fid, '%e ',total_u(1,end));fprintf(fid, '],');
		            fprintf(fid,'\n   ['); fprintf(fid, '%e, ',total_u(2,2:end-1)); fprintf(fid, '%e ',total_u(2,end));fprintf(fid, '],');
		            fprintf(fid,'\n   ['); fprintf(fid, '%e, ',total_u(3,2:end-1)); fprintf(fid, '%e ',total_u(3,end));fprintf(fid, '],');
		            fprintf(fid,'\n   ['); fprintf(fid, '%e, ',total_u(1,1:end-2)); fprintf(fid, '%e ',total_u(1,end-1));fprintf(fid, '],');
		            fprintf(fid,'\n   ['); fprintf(fid, '%e, ',total_u(2,1:end-2)); fprintf(fid, '%e ',total_u(2,end-1));fprintf(fid, '],');
		            fprintf(fid,'\n   ['); fprintf(fid, '%e, ',total_u(3,1:end-2)); fprintf(fid, '%e ',total_u(3,end-1));fprintf(fid, ']\n');
		            fprintf(fid, '  ]\n');
		            fprintf(fid, '}\n');
		            fclose(fid);
		            // fprintf('(%d) total_loading: %e\n',sample_id,total_u(end));
		            
		            // Time step
		            fid = fopen([dir_name,'/timesteps.json'], 'w');
		            fprintf(fid, '{\n');
		            fprintf(fid, '  "number_of_time_steps": ');
		            fprintf(fid,'%d, ',n_micro_step);
		            fprintf(fid, '\n  "time_steps": [ 0.0, ');
		            fprintf(fid,'%e, ', total_t(1:end-1));
		            fprintf(fid,'%e',total_t(end));
		            fprintf(fid, '],\n  "number_of_steps_to_print": ');
		            fprintf(fid,'%d, ',numel(print_time_step));
		            fprintf(fid, '\n  "list_of_steps_to_print": [');
		            
		            switch opt_ts_print
		                case 0 // all steps
		                    fprintf(fid,'%d, ',print_time_step(1:end-1));
		                    fprintf(fid,'%d',print_time_step(end));
		                case 1 // the first & last steps
		                    fprintf(fid,'%d, %d', print_time_step(1), print_time_step(end));
		                case 2
		                    fprintf(fid,'%d, ',print_time_step(1:end-1));
		                    fprintf(fid,'%d',print_time_step(end));
		                otherwise 
		                    fprintf("ERROR: please specify option for the number of time steps printed\r\n");
		                    return
		            end
		            
		            fprintf(fid, ']\n}\n');
		            fclose(fid);
		            fprintf('(%d) total_time: %e\n',sample_id,total_t(end));
		            
		            // Mechanical bcs
		            fid = fopen([dir_name,'/Mechanical_bc.json'], 'w');
		            fprintf(fid,'{\n');
		            fprintf(fid,'  "#": "boundary condition description: geom_type, geom_id, flag(u), flag(v), flag(w)",\n');
		            fprintf(fid,'  "#": "first 2 will always be there. Additional elements per line are based on the # of degrees of freedom in multiphysics.in",\n');
		            fprintf(fid,'  "#": "replace -1 with this value",\n');
		            fprintf(fid,'  "bc_data":    [\n');
		            fprintf(fid,'                [1, 69, 1, 1, 1 ],\n');
		            fprintf(fid,'                [1, 70, 1, 1, 1 ],\n');
		            fprintf(fid,'                [1, 71, 1, 1, 1 ],\n');
		            fprintf(fid,'                [1, 72, 1, 1, 1 ],\n');
		            fprintf(fid,'                [1, 73, 1, 1, 1 ],\n');
		            fprintf(fid,'                [1, 74, 1, 1, 1 ],\n');
		            fprintf(fid,'                [1, 75, 1, 1, 1 ],\n');
		            fprintf(fid,'                [1, 76, 1, 1, 1 ],\n');
		            fprintf(fid,'                [2, 125, 1, 1, 1 ],\n');
		            fprintf(fid,'                [2, 126, 1, 1, 1 ],\n');
		            fprintf(fid,'                [2, 127, 1, 1, 1 ],\n');
		            fprintf(fid,'                [2, 128, 1, 1, 1 ],\n');
		            fprintf(fid,'                [2, 129, 1, 1, 1 ],\n');
		            fprintf(fid,'                [2, 130, 1, 1, 1 ],\n');
		            fprintf(fid,'                [2, 131, 1, 1, 1 ],\n');
		            fprintf(fid,'                [2, 132, 1, 1, 1 ],\n');
		            fprintf(fid,'                [5, 23, 1, 1, 1 ],\n');
		            fprintf(fid,'                [5, 26, 1, 1, 1 ]\n');
		            fprintf(fid,'                ],\n');
		            fprintf(fid,'    "#": "replace -1 with this value",\n');
		            fprintf(fid,'    "replacements": [ ');
		            fprintf(fid,'%3.5e, %3.5e, %3.5e, 0.0, 0.0, 0.0',d_jump(1),d_jump(2),d_jump(3));
		            fprintf(fid,' ]\n');
		            fprintf(fid,'}\n');
		            fclose(fid);                     
		        case 1 // parallelized input type
		            // STEP 3: parallelized inputs, utilizing in.st
		            fid = fopen([dir_name,'/Mechanical.bcv'], 'w');
		            fprintf(fid,'6\n');
		            fprintf(fid,'%3.5e\n%3.5e\n%3.5e\n0.0\n0.0\n0.0\n',d_jump(1),d_jump(2),d_jump(3));
		            fclose(fid);
		            
		            // NOTE: MATLAB read the column first
		            load_matrix = [total_u(1,2:end);
		                           total_u(2,2:end);
		                           total_u(3,2:end);
		                           total_u(1,1:end-1);
		                           total_u(2,1:end-1);
		                           total_u(3,1:end-1)];
		                       
		            fid = fopen([dir_name,'/Mechanical.load'], 'w');
		            fprintf(fid,'%d\n',n_micro_step-1);
		            fprintf(fid,'%d ',1:1:n_micro_step-1);
		            fprintf(fid,'\n\n');
		            fprintf(fid,'%3.5e %3.5e %3.5e %3.5e %3.5e %3.5e ',load_matrix);
		            fclose(fid);
		            
		            
		            fid = fopen([dir_name,'/in_st_p/',file_base,'_0.in.st'], 'w');
		            fprintf(fid,'1e-05 10 4 1\n');
		            fprintf(fid,'\n');
		            fprintf(fid,'%d\n',n_micro_step);
		            fprintf(fid,'0.0 ');
		            fprintf(fid,'%f ',total_t(1:end));
		            fprintf(fid,'\n\n');
		            fprintf(fid,'%d\n',numel(print_time_step));
		            switch opt_ts_print
		                case 0 // all steps
		                    fprintf(fid,'%d ',print_time_step(1:end));
		                case 1 // the first & last steps
		                    fprintf(fid,'%d %d', print_time_step(1), print_time_step(end));
		                case 2
		                    fprintf(fid,'%d ',print_time_step(1:end));
		                otherwise
		                    fprintf("ERROR: please specify option for the number of time steps printed\r\n");
		                    return
		            end
		            fclose(fid);
		            
		            
		            fid = fopen([dir_name,'/in_st_t/',file_base,'_0.in.st'], 'w');
		            fprintf(fid,'1e-05 10 4 5\n');
		            fprintf(fid,'\n');
		            fprintf(fid,'%d\n',n_micro_step);
		            fprintf(fid,'0.0 ');
		            fprintf(fid,'%f ',total_t(1:end));
		            fprintf(fid,'\n\n');
		            fprintf(fid,'%d\n',numel(print_time_step));            
		            
		            switch (opt_ts_print) {
		                case 0 // all steps
		                    fprintf(fid,'%d ',print_time_step(1:end));
		                    break;
		                case 1 // the first & last steps
		                    fprintf(fid,'%d %d', print_time_step(1), print_time_step(end));
		                    break;
		                case 2
		                    fprintf(fid,'%d ',print_time_step(1:end));                    
		                    break;
		                default:
		                    fprintf("ERROR: please specify option for the number of time steps printed\r\n");
		                    break;
		            }   
		            fclose(fid);            
		        default:
		            fprintf("ERROR: please specify desired inputs type\n");
		            break;
		    end // end of switch opt_par
		} // end of sample_id 

	return err;

}