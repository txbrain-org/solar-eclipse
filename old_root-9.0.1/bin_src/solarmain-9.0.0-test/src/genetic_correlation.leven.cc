#include <Eigen/Dense>
#include <iostream>
#include "solar.h"
#include <fstream>
#include "solar_mle_setup.h"
#include <cmath>
#include <string>
#include <iomanip>
#include <chrono>
#include <stdio.h>
#include <cstdlib>
#include <omp.h>
using namespace std;
#define MAX_ITERATIONS 500
#define MAX_DELTA_ERROR 1e-07
#define MAX_LOGLIK_ERROR 1e-09
#define DP 1.0E-7
#define STOP_CRITERIA 1.0E-3
#define T 1.0

extern "C" void cdfchi_ (int*, double*, double*, double*, double*,
                         int*, double*);
static double chicdf(double chi, double df){
    double p, q, bound;
    int status = 0;
    int which = 1;
    
    
    cdfchi_ (&which, &p, &q, &chi, &df, &status, &bound);
    
    return q/2.0;
}
static inline double calculate_constraint(const double x){
    //return exp(x)/(1 + exp(x));
    return x*x/(1.0 + x*x);

}
static double reverse_constraint(double x){
    return sqrt(x/(1-x));
    
    
}
static inline double calculate_dconstraint(const double x){

    return 2*x/pow(1+x*x, 2);
   // const double e_x = exp(x);
	//return e_x*pow(e_x + 1, -2);
	
}


static inline double calculate_ddconstraint(const double x){

    return -2*(3*x*x - 1)/pow((x*x + 1), 3);
	//const double e_x = exp(x);
	//return -(e_x-1.0)*e_x*pow(e_x+1,-3);
	
}
static void compute_blockwise_diagonal_inversions(Eigen::VectorXd & A, Eigen::VectorXd & BC, \
                                                     Eigen::VectorXd & D){
/*
   Eigen::VectorXd denom = A.cwiseProduct(D) - BC.cwiseAbs2();
    denom = denom.cwiseInverse();
    BC = -BC.cwiseProduct(denom);
    Eigen::VectorXd old_A = A;
    A = D.cwiseProduct(denom);
    D = old_A.cwiseProduct(denom);*/
    Eigen::VectorXd A_inverse = A.cwiseInverse();
    Eigen::VectorXd D_inverse = D.cwiseInverse();
    
    Eigen::VectorXd E = (D-BC.cwiseAbs2().cwiseProduct(A_inverse)).cwiseInverse();
        
    Eigen::VectorXd new_A = A_inverse + A_inverse.cwiseAbs2().cwiseProduct(BC.cwiseAbs2().cwiseProduct(E));
    Eigen::VectorXd new_BC = -A_inverse.cwiseProduct(BC.cwiseProduct(E));
    Eigen::VectorXd new_D = E;
    A = new_A;
    BC = new_BC;
    D = new_D;
}
static double calculate_rho(const double x){
   // return x;
    return atan(x)*2.0/M_PI;
}
static double reverse_rho(double x){
        return tan(x)*2.0/M_PI;
}
static double calculate_rho_dconstraint(const double x){
    return 2.0/((1.0 + x*x)*M_PI);
   
}
static double calculate_rho_ddconstraint(const double x){
    return -4.0*x/(pow((1.0 + x*x), 2)*M_PI);
   
}
static Eigen::VectorXd combine_parameters_and_beta(Eigen::VectorXd parameters, Eigen::VectorXd beta){
	Eigen::VectorXd all_parameters(parameters.rows() + beta.rows());
	for(int i = 0; i < parameters.rows(); i++){
		all_parameters(i) = parameters(i);
	}

	for(int i = 0; i < beta.rows(); i++){
		all_parameters(parameters.rows() + i) = beta(i);
	}
	return all_parameters;
}
static double calculate_loglikelihood_param(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two, Eigen::MatrixXd covariate_matrix, Eigen::VectorXd lambda, Eigen::VectorXd parameters, const bool use_constraints = true){
                   
    Eigen::VectorXd omega_one_one;// = omega_1_1 = sd_one*sd_one*(h2r_one*lambda + (1.0-h2r_one)*Eigen::VectorXd::Ones(lambda.rows()));
    Eigen::VectorXd omega_one_two;
     Eigen::VectorXd omega_two_two;
     const double h2r_one = (use_constraints) ? calculate_constraint(parameters(0)) : parameters(0);
     const double h2r_two = (use_constraints) ? calculate_constraint(parameters(1)) : parameters(1);
     const double rhog = (use_constraints) ? calculate_rho(parameters(4)) : parameters(4);
     const  double rhoe =  (use_constraints) ? calculate_rho(parameters(5)) : parameters(5);
     const double sd_one = fabs(parameters(2));
     const  double sd_two = fabs(parameters(3));
      omega_one_one =  sd_one*sd_one*(h2r_one*lambda + (1.0-h2r_one)*Eigen::VectorXd::Ones(lambda.rows()));
      omega_one_two = sd_one*sd_two*(lambda*sqrt(h2r_one)*sqrt(h2r_two)*rhog + rhoe*sqrt(1.0-h2r_one)*sqrt(1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
      omega_two_two = sd_two*sd_two*(h2r_two*lambda + (1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
    

    Eigen::VectorXd beta_one(covariate_matrix.cols());
    Eigen::VectorXd beta_two(covariate_matrix.cols());
    for(int i = 0 ; i < covariate_matrix.cols(); i++){
        beta_one(i) = parameters(i + 6);
        beta_two(i) = parameters(i + covariate_matrix.cols() + 6);
    }
    Eigen::VectorXd residual_one = Y_one - covariate_matrix*beta_one;//parameters(0);
    Eigen::VectorXd residual_two = Y_two - covariate_matrix*beta_two;//parameters(1);

     
    
    double part_one = 0.0;

   Eigen::VectorXd omega_det = omega_one_one.cwiseProduct(omega_two_two) - omega_one_two.cwiseAbs2();
    for(int i = 0; i < lambda.rows() ; i++){
        part_one += log(fabs(omega_det(i)));
    }

    compute_blockwise_diagonal_inversions(omega_one_one, omega_one_two, omega_two_two);
   

    const double part_two = residual_one.dot(residual_one.cwiseProduct(omega_one_one)) + residual_two.dot(residual_two.cwiseProduct(omega_two_two)) + 2.0*residual_one.dot(omega_one_two.cwiseProduct(residual_two));
    if(use_constraints)
        

    return -0.5*(part_one + part_two);
}

static double calculate_gradient(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two,\
                        Eigen::MatrixXd covariate_matrix, Eigen::VectorXd lambda, \
                        Eigen::VectorXd all_parameters, const int index,  bool is_hessian = true){

    double h;
   // is_hessian = false;
    if(is_hessian){
        h = pow(DP, 0.66667)*max(abs(all_parameters(index)), 1.0);
    }else{
        h = DP*max(abs(all_parameters(index)), 1.0);
    }
    Eigen::VectorXd positive_parameters = all_parameters;
    positive_parameters(index) += h;
    Eigen::VectorXd negative_parameters =all_parameters;
    negative_parameters(index) -= h;
    return (calculate_loglikelihood_param(Y_one,  Y_two,  covariate_matrix, \
     lambda,  positive_parameters) - calculate_loglikelihood_param(Y_one,  \
     Y_two,  covariate_matrix, lambda,  negative_parameters))/(2.0*h);
}
static double calculate_hessian(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two,\
                        Eigen::MatrixXd covariate_matrix, Eigen::VectorXd lambda, \
                        Eigen::VectorXd all_parameters, const int index_one,\
                        const int index_two){
    Eigen::VectorXd positive_parameters = all_parameters;
    double h = pow(DP, 1)*max(abs(all_parameters(index_two)), 1.0);
    positive_parameters(index_two) += h;
    
    return ((calculate_gradient(Y_one, Y_two, covariate_matrix, \
        lambda, positive_parameters,  index_one, false)) - calculate_gradient(Y_one, Y_two, covariate_matrix, \
        lambda, all_parameters,  index_one, false))/(h);
}

static void genetic_correlation_calculate_hessian_and_gradient(Eigen::MatrixXd & hessian, Eigen::VectorXd & gradient, Eigen::VectorXd Y_one, Eigen::VectorXd Y_two,\
                                                                Eigen::MatrixXd covariate_matrix,  Eigen::VectorXd lambda, Eigen::VectorXd all_parameters){
                                                                




  	for(int index = 0; index < all_parameters.rows() ; index++){

        gradient(index) = calculate_gradient(Y_one, Y_two,\
                        covariate_matrix, lambda, \
                        all_parameters, index);

        hessian(index, index) = calculate_hessian(Y_one, Y_two,\
                        covariate_matrix, lambda, \
                        all_parameters, index, \
                         index);

    }
  

    for(int i = 0 ; i < all_parameters.rows(); i++){

        for(int j = i  + 1; j < all_parameters.rows() ; j++){

         hessian(i, j) = hessian(j, i) = calculate_hessian(Y_one, Y_two,\
                        covariate_matrix, lambda, \
                        all_parameters, i, j);
        }
       
    }                                                                  
                                                                
} 





static void calculate_mean_and_sd(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two, Eigen::MatrixXd covariate_matrix, Eigen::VectorXd lambda, Eigen::VectorXd & parameters, Eigen::VectorXd & beta){
    Eigen::VectorXd omega_one = lambda*calculate_constraint(parameters(0)) + (1.0 - calculate_constraint(parameters(0)))*Eigen::VectorXd::Ones(Y_one.rows());
    Eigen::VectorXd omega_two = lambda*calculate_constraint(parameters(1)) + (1.0 - calculate_constraint(parameters(1)))*Eigen::VectorXd::Ones(Y_one.rows());

    omega_one = omega_one.cwiseInverse();
    
    Eigen::VectorXd beta_one = (covariate_matrix.transpose()*omega_one.asDiagonal()*covariate_matrix).inverse()*covariate_matrix.transpose()*omega_one.asDiagonal()*Y_one;

    parameters(2) = sqrt((Y_one - covariate_matrix*beta_one).cwiseAbs2().dot(omega_one)/Y_one.rows());
    
    omega_two = omega_two.cwiseInverse();

    Eigen::VectorXd beta_two = (covariate_matrix.transpose()*omega_two.asDiagonal()*covariate_matrix).inverse()*covariate_matrix.transpose()*omega_two.asDiagonal()*Y_two;
 
  
    parameters(3) = sqrt((Y_two - covariate_matrix*beta_two).cwiseAbs2().dot(omega_two)/Y_two.rows());
    for(int i = 0; i < covariate_matrix.cols(); i++){
        beta(i) = beta_one(i);
        beta(i + covariate_matrix.cols()) = beta_two(i);
    }    

}

static double calculate_bivariate_model(Eigen::VectorXd trait_one, Eigen::VectorXd trait_two, Eigen::MatrixXd covariate_matrix, \
    Eigen::VectorXd eigenvalues, Eigen::VectorXd & final_parameters, Eigen::VectorXd & final_beta, Eigen::VectorXd & standard_errors, \
    double  & best_h, int & iteration_count, int constrain_parameter = 0, bool debug = false){

    double best_loglikelihood;
    bool first_loglik = false;
    double h;
    Eigen::VectorXd beta(covariate_matrix.cols()*2);
    Eigen::VectorXd parameters(6);
    Eigen::VectorXd all_parameters(6 + beta.rows());
    parameters(0) = 1.0;
    parameters(1) = 1.0;
    parameters(4) = 0.0;
    parameters(5) = 0.0;

    calculate_mean_and_sd(trait_one, trait_two, covariate_matrix, eigenvalues, parameters ,beta);
    all_parameters = combine_parameters_and_beta(parameters, beta);
    double loglik = calculate_loglikelihood_param(trait_one, trait_two, covariate_matrix, eigenvalues, all_parameters);

    Eigen::VectorXd delta = Eigen::VectorXd::Zero(6 + covariate_matrix.cols()*2);
    Eigen::MatrixXd hessian(6 + covariate_matrix.cols()*2, 6 + covariate_matrix.cols()*2);
    Eigen::VectorXd gradient(6 + covariate_matrix.cols()*2);
    int iteration = 0;
    double last_loglik = 0.0;
    double loglik_error;
    double max_delta;
    genetic_correlation_calculate_hessian_and_gradient(hessian,  gradient, trait_one, trait_two,\
                                                               covariate_matrix,  eigenvalues,  all_parameters);


    if(constrain_parameter == 1){
        gradient(4) = 0.0;
        for(int i = 0 ;i < hessian.rows(); i++){
            hessian(i, 4) = hessian(4, i) = 0.0;
        }
        hessian(4, 4) = 1.0;
    }
    if(constrain_parameter == 2){
        gradient(5) = 0.0;
        for(int i = 0 ;i < hessian.rows(); i++){
            hessian(i, 5) = hessian(5, i) = 0.0;
        }
        hessian(5, 5) = 1.0;
    }
    /*
    if(constrain_parameter == 1){
    	Eigen::MatrixXd new_hessian(5 + beta.rows(), 5 + beta.rows());
    	Eigen::VectorXd new_gradient(5 + beta.rows());
   

    	for(int i = 0 ; i < hessian.rows(); i++){
    		int index_i = i;
    		if(i == 4){
    			continue;
    		}else if (i > 4){
    			index_i = i -1;
    		}
    		new_gradient(index_i) = gradient(i);
    		for(int j = 0 ; j < hessian.rows(); j++){
    			int index_j = j;
    			if(j == 4){
    				continue;
    			}else if (j > 4){
    				index_j = j -1;
    			}
    			new_hessian(index_i, index_j) = hessian(i, j);	
    		}
    	}
    	Eigen::VectorXd new_delta = -new_hessian.inverse()*new_gradient;
    	for(int i = 0 ; i < delta.rows(); i++){
    		int index_i = i;
    		if(i == 4){
    			delta(i) = 0.0;
    			continue;
    		}else if (i > 4){
    			index_i = i - 1;
    		}
    		delta(i) = new_delta(index_i);
    	}	
    }else if(constrain_parameter == 2){
    	Eigen::MatrixXd new_hessian(5 + beta.rows(), 5 + beta.rows());
    	Eigen::VectorXd new_gradient(5 + beta.rows());
    	

    	for(int i = 0 ; i < hessian.rows(); i++){
    		int index_i = i;
    		if(i == 5){
    			continue;
    		}else if (i > 5){
    			index_i = i -1;
    		}
    		new_gradient(index_i) = gradient(i);
    		for(int j = 0 ; j < hessian.rows(); j++){
    			int index_j = j;
    			if(j == 5){
    				continue;
    			}else if (j > 5){
    				index_j = j -1;
    			}
    			new_hessian(index_i, index_j) = hessian(i, j);	
    		}
    	}
    	Eigen::VectorXd new_delta = -new_hessian.inverse()*new_gradient;
    	for(int i = 0 ; i < delta.rows(); i++){
    		int index_i = i;
    		if(i == 5){
    			delta(i) = 0.0;
    			continue;
    		}else if (i > 5){
    			index_i = i - 1;
    		}
    		delta(i) = new_delta(index_i);
    	}
    }else{                                                                           	

        delta = -hessian.inverse()*gradient;
    }*/

   // delta = -hessian.inverse()*gradient;
    delta = -hessian.inverse()*gradient;
    double error;
    int converge_count = 0;
    do{
        
    	Eigen::VectorXd test_parameters = all_parameters + delta;
    	//test_parameters += delta;
    	double best_loglik = calculate_loglikelihood_param(trait_one, trait_two, covariate_matrix, eigenvalues, test_parameters);
    	Eigen::VectorXd best_parameters = test_parameters;
    	for(int step = 1; step <= 10; step++){
            Eigen::MatrixXd new_hessian = hessian + step*0.1*hessian.diagonal().asDiagonal();
            if(constrain_parameter == 1) new_hessian(4, 4) = 1.0;
            if(constrain_parameter == 2) new_hessian(5, 5) = 1.0;
            Eigen::VectorXd test_delta = -new_hessian*gradient;
    		test_parameters = all_parameters + test_delta;
    		double test_loglik = calculate_loglikelihood_param(trait_one, trait_two, covariate_matrix, eigenvalues, test_parameters);
    		if(test_loglik > best_loglik || (best_loglik != best_loglik && test_loglik == test_loglik)){
    			best_loglik = test_loglik;
    			best_parameters = test_parameters;
    		}
    	}
        last_loglik = loglik;
        loglik = best_loglik;
        loglik_error = fabs((last_loglik - loglik)/last_loglik);
        all_parameters = best_parameters;
    	genetic_correlation_calculate_hessian_and_gradient(hessian,  gradient, trait_one, trait_two,\
                                                               covariate_matrix,  eigenvalues,  all_parameters);
        if(constrain_parameter == 1){
            gradient(4) = 0.0;
            for(int i = 0 ;i < hessian.rows(); i++){
                hessian(i, 4) = hessian(4, i) = 0.0;
            }
            hessian(4, 4) = 1.0;
        }
        if(constrain_parameter == 2){
            gradient(5) = 0.0;
            for(int i = 0 ;i < hessian.rows(); i++){
                hessian(i, 5) = hessian(5, i) = 0.0;
            }
            hessian(5, 5) = 1.0;
        }/*    
    	if(constrain_parameter == 1){
    		Eigen::MatrixXd new_hessian(5 + beta.rows(), 5 + beta.rows());
    		Eigen::VectorXd new_gradient(5 + beta.rows());
   

    		for(int i = 0 ; i < hessian.rows(); i++){
    			int index_i = i;
    			if(i == 4){
    				continue;
    			}else if (i > 4){
    				index_i = i -1;
    			}
    			new_gradient(index_i) = gradient(i);
    			for(int j = 0 ; j < hessian.rows(); j++){
    			    int index_j = j;
    				if(j == 4){
    					continue;
    				}else if (j > 4){
    					index_j = j -1;
    				}
    				new_hessian(index_i, index_j) = hessian(i, j);	
    			}
    		}
    		Eigen::VectorXd new_delta = -new_hessian.inverse()*new_gradient;
    		for(int i = 0 ; i < delta.rows(); i++){
    			int index_i = i;
    			if(i == 4){
    				delta(i) = 0.0;
    				continue;
    			}else if (i > 4){
    				index_i = i - 1;
    			}
    			delta(i) = new_delta(index_i);
    		}	
    	}else if(constrain_parameter == 2){
    		Eigen::MatrixXd new_hessian(5 + beta.rows(), 5 + beta.rows());
    		Eigen::VectorXd new_gradient(5 + beta.rows());
    	

    		for(int i = 0 ; i < hessian.rows(); i++){
    			int index_i = i;
    			if(i == 5){
    				continue;
    			}else if (i > 5){
    				index_i = i -1;
    			}
    			new_gradient(index_i) = gradient(i);
    			for(int j = 0 ; j < hessian.rows(); j++){
    			    int index_j = j;
    				if(j == 5){
    					continue;
    				}else if (j > 5){
    					index_j = j -1;
    				}
    				new_hessian(index_i, index_j) = hessian(i, j);	
    			}
    		}
    		Eigen::VectorXd new_delta = -new_hessian.inverse()*new_gradient;
    		for(int i = 0 ; i < delta.rows(); i++){
    			int index_i = i;
    			if(i == 5){
    				delta(i) = 0.0;
    				continue;
    			}else if (i > 5){
    				index_i = i - 1;
    			}
    			delta(i) = new_delta(index_i);
    		}
    	}else{                                                                           	

        	delta = -hessian.inverse()*gradient;
        }*/
        delta = -hessian.inverse()*gradient;
       // error = fabs(delta.dot(gradient));
        if (loglik_error <= MAX_LOGLIK_ERROR){
            converge_count++;
        }else{
            converge_count = 0;
        } 
    
        iteration++; 
    }while(loglik == loglik && delta == delta && iteration < MAX_ITERATIONS  && (converge_count  < 3));  
    if(loglik == loglik && delta == delta && (loglik > best_loglikelihood || first_loglik == false)){

        
        best_loglikelihood = loglik;
        for(int i = 0 ; i < 6 ; i++){
        	parameters(i) = all_parameters(i);
        }
        for(int i = 0 ; i < beta.rows(); i++){
        	beta(i) = all_parameters(6 + i);
        }
        final_parameters = parameters;
        final_beta = beta;
        first_loglik = true;
        iteration_count = iteration;

    	standard_errors = (-hessian).inverse().diagonal().cwiseAbs().cwiseSqrt();
    	standard_errors(0) =fabs(standard_errors(0)*calculate_dconstraint(all_parameters(0)));
    	standard_errors(1) =fabs(standard_errors(1)*calculate_dconstraint(all_parameters(1)));
    	standard_errors(4) =(constrain_parameter == 1) ? 0.0 : fabs(standard_errors(4)*calculate_rho_dconstraint(all_parameters(4)));
    	standard_errors(5) =(constrain_parameter == 2) ? 0.0 : fabs(standard_errors(5)*calculate_rho_dconstraint(all_parameters(5)));          
      
    }
    
  
  	if(first_loglik == false){
    	standard_errors = Eigen::VectorXd::Zero(6+covariate_matrix.cols()*2);
    	return nan("");
  	}

  
	return best_loglikelihood;
}                                                             
 
static const char * calculate_genetic_correlation(Tcl_Interp * interp, vector<string> trait_list, const char* phenotype_filename, bool display_pvalues, bool debug, double main_best_h, const char * evd_data_filename = 0){
    vector<string> cov_list;
    vector<string> unique_cov_terms;
    int success;
    Covariate * c;
    
    int n_covariates = 0;
    for (int i = 0;( c = Covariate::index(i)); i++)
    {
        char buff[512];
        c->fullname(&buff[0]);
        cov_list.push_back(string(&buff[0]));
        CovariateTerm * cov_term;
        
        for(cov_term = c->terms(); cov_term; cov_term = cov_term->next){
            bool found = false;
            
            for(vector<string>::iterator cov_iter = unique_cov_terms.begin(); cov_iter != unique_cov_terms.end(); cov_iter++){
                if(!StringCmp(cov_term->name, cov_iter->c_str(), case_ins)){
                    found = true;
                    break;
                }
            }
            if(!found){
                unique_cov_terms.push_back(string(cov_term->name));
            }
        }
        n_covariates++;
    }
    vector<string> field_list;
    field_list.push_back(trait_list[0]);
    field_list.push_back(trait_list[1]);
    for(int i = 0; i < unique_cov_terms.size(); i++){
        field_list.push_back(unique_cov_terms[i]);
    }    
    solar_mle_setup * file_data;
    try{
        file_data = new solar_mle_setup(field_list, phenotype_filename,interp,  true); 
    }catch(Parse_Expression_Error & error){
    	
      
        return error.what().c_str();

    }catch(Solar_File_Error & error){
        return error.what().c_str();
 
    }catch(Expression_Eval_Error & error){
        return error.what().c_str();
      
    }catch(Misc_Error & error){
        
        
        return error.what().c_str();
    }catch(Syntax_Error &e){
    
    	return "Syntax Error occurred in expression";

    }catch(Undefined_Function &e){
    	return "Undefined Function Error occurred in expression";

    }catch(Unresolved_Name &e){
    	return "Unresolved Name Error occurred in expression";
 
    }catch(Undefined_Name &e){
        cout << e.name << endl;
    	return "Undefined Name Error occurred in expression";
    }catch(...){
    	
    	return "Unkown error occurred reading phenotype data";
    }
    Eigen::MatrixXd field_matrix = file_data->return_output_matrix();
    Eigen::VectorXd trait_one =  field_matrix.col(0);
    Eigen::VectorXd trait_two =  field_matrix.col(1);
    Eigen::MatrixXd covariate_matrix = Eigen::ArrayXXd::Ones(trait_one.rows(), n_covariates + 1);
    if (n_covariates > 0){
        Eigen::MatrixXd covariate_term_matrix(trait_one.rows(), field_list.size() - 2);
   
        for(int col = 0; col < unique_cov_terms.size(); col++){
        
            if(!StringCmp(unique_cov_terms[col].c_str(), "SEX", case_ins)){
                if((field_matrix.col(col + 2).array() == 2.0).count() != 0){
                for(int row = 0 ; row < covariate_term_matrix.rows(); row++){
                    if(field_matrix(row, col + 2) == 2.0){
                        covariate_term_matrix(row, col) = 1.0;
                    }else{
                        covariate_term_matrix(row, col) = 0.0;
                    }
                }
            }else{
                covariate_term_matrix.col(col) = field_matrix.col(col + 2);
            }
            
            }else if(strstr(unique_cov_terms[col].c_str(), "snp_") != NULL || strstr(unique_cov_terms[col].c_str(), "SNP_") != NULL){
                covariate_term_matrix.col(col) = field_matrix.col(col + 2).array() - field_matrix.col(col + 2).mean();
                continue;
            } else {
            covariate_term_matrix.col(col) =  field_matrix.col(col + 2).array() - field_matrix.col(col + 2).mean();
            }
        }
    
        Covariate * cov;
    
        for(int col = 0; (cov = Covariate::index(col)); col++){
            CovariateTerm * cov_term;
            for(cov_term = cov->terms(); cov_term; cov_term = cov_term->next){
                int index = 0;
            
                for(vector<string>::iterator cov_iter = unique_cov_terms.begin(); cov_iter != unique_cov_terms.end(); cov_iter++){
                    if(!StringCmp(cov_term->name, cov_iter->c_str(), case_ins)){
                        break;
                    }
                    index++;
                }
            
                if(cov_term->exponent == 1){
                    covariate_matrix.col(col) = covariate_matrix.col(col).array()*covariate_term_matrix.col(index).array();
                }else{
                    covariate_matrix.col(col) = covariate_matrix.col(col).array()*pow(covariate_term_matrix.col(index).array(), cov_term->exponent);
                }
            }
        
        }
    }    
    Eigen::MatrixXd eigenvectors; 
    Eigen::VectorXd eigenvalues;
    if(evd_data_filename){
        vector<string> evd_ids;
        string evd_id_filename = string(evd_data_filename) + string(".ids");
        ifstream ids_in(evd_id_filename.c_str());
       if(ids_in.is_open() == false){
            string error_string = "File " + evd_id_filename + ".ids not found";
            delete file_data;
            return error_string.c_str();
        }         
        string current_id;
        while (ids_in >> current_id) {
            evd_ids.push_back(current_id);
        }
        ids_in.close();
        vector<string> mle_setup_ids = file_data->get_ids();
        if(mle_setup_ids.size() != evd_ids.size()){
            string error_string = "ID count between EVD data IDs and FPHI trait reader IDs doesn't match\n EVD data ID count: " + to_string(evd_ids.size()) + " FPHI trait reader ID count: " + to_string(mle_setup_ids.size());
            delete file_data;
            return error_string.c_str();
        }
        for(int i  = 0 ; i < evd_ids.size(); i++){
            if(evd_ids[i] != mle_setup_ids[i]){
                string error_string = "At ID index " + to_string(i) + " EVD ID " + evd_ids[i] + "does not match FPHI ID " + mle_setup_ids[i];
            }
       }
       eigenvectors = Eigen::MatrixXd::Zero(evd_ids.size(), evd_ids.size());
       eigenvalues = Eigen::VectorXd::Zero(evd_ids.size());
       string evd_eigenvectors_filename = string(evd_data_filename) + ".eigenvectors"; 
       string evd_eigenvalues_filename = string(evd_data_filename) + ".eigenvalues";
       ifstream eigenvectors_stream(evd_eigenvectors_filename.c_str());
       if(eigenvectors_stream.is_open() == false){
            string error_string = "File " + evd_eigenvectors_filename  + ".eigenvectors not found";
            delete file_data;
            return error_string.c_str();
        }
       ifstream eigenvalues_stream(evd_eigenvalues_filename.c_str());
       if(eigenvalues_stream.is_open() == false){
            string error_string = "File " + evd_eigenvalues_filename + ".eigenvalues not found";
            delete file_data;
            return error_string.c_str();
        }       
       for(int row = 0 ; row < evd_ids.size(); row++){
            eigenvalues_stream >> eigenvalues(row) ;
            for(int col = 0; col < evd_ids.size() ; col++){
                eigenvectors_stream >>  eigenvectors(row, col);
            }
      }
    
      eigenvectors_stream.close();
      eigenvalues_stream.close();
    }else{
                 
        eigenvectors = file_data->get_eigenvectors().transpose();
        eigenvalues = file_data->get_eigenvalues();
    }
    Eigen::VectorXd ones = Eigen::ArrayXd::Ones(trait_one.rows());
    Eigen::VectorXd S_ones = eigenvectors*ones;
    cout.precision(7);
    trait_one = eigenvectors*trait_one;
    trait_two = eigenvectors*trait_two;
    covariate_matrix = eigenvectors*covariate_matrix;


    Eigen::VectorXd final_parameters(6);
    Eigen::VectorXd final_beta(covariate_matrix.cols()*2);
   // double main_best_h;
    int main_iteration_count;
    Eigen::VectorXd standard_errors;
    double main_loglik = calculate_bivariate_model(trait_one,  trait_two,  covariate_matrix, eigenvalues,  final_parameters, final_beta, standard_errors, main_best_h, main_iteration_count, 0,  debug); 
    if(main_loglik != main_loglik ){
        Solar_Eval(interp, "loglike set 0");
        const char * error = "Convergence of this model could not be achieved";
        delete file_data;
        return error;
    }
    string loglik_command = "loglike set " + to_string(main_loglik);
    Solar_Eval(interp, loglik_command.c_str());   
    if(display_pvalues){
        Eigen::VectorXd rhog_parameters(6);
        Eigen::VectorXd rhog_beta(covariate_matrix.cols()*2);
        int iteration_count;
        double rhog_best_h = main_best_h;
        Eigen::VectorXd rhog_standard_errors;
        double rhog_loglik  = calculate_bivariate_model(trait_one,  trait_two,  covariate_matrix, eigenvalues,  rhog_parameters, rhog_beta, rhog_standard_errors, rhog_best_h, main_iteration_count, 1,  debug); 
        
        
        
        Eigen::VectorXd rhoe_parameters(6);
        Eigen::VectorXd rhoe_beta(covariate_matrix.cols()*2);
        double rhoe_best_h = main_best_h;
        Eigen::VectorXd rhoe_standard_errors;
        double rhoe_loglik  = calculate_bivariate_model(trait_one,  trait_two,  covariate_matrix, eigenvalues,  rhoe_parameters, rhoe_beta, rhoe_standard_errors, rhoe_best_h, main_iteration_count, 2,  debug); 

       
        double rhog_pvalue = 1.0;
        if(rhog_loglik == rhog_loglik){
            rhog_pvalue = 2.0*chicdf(2.0*(main_loglik - rhog_loglik), 1); 
        }
        double rhoe_pvalue = 1.0;
        if(rhoe_loglik == rhoe_loglik){
            rhoe_pvalue = 2.0*chicdf(2.0*(main_loglik - rhoe_loglik), 1);
        }
        const double rhop = calculate_rho(final_parameters(4))*sqrt(calculate_constraint(final_parameters(0)))*\
             sqrt(calculate_constraint(final_parameters(1))) + calculate_rho(final_parameters(5))*\
             sqrt(1.0-calculate_constraint(final_parameters(0)))*sqrt(1.0-calculate_constraint(final_parameters(1)));
        final_parameters(0) = calculate_constraint(final_parameters(0));
        final_parameters(1) = calculate_constraint(final_parameters(1)); 
        final_parameters(4) = calculate_rho(final_parameters(4));
        final_parameters(5) = calculate_rho(final_parameters(5));               
        cout << endl;     
        //string solar_command = "set RHOP " + to_string(rhop);
        cout <<  "*******************************" << endl;
        cout <<  "*     Genetic Correlation     *" << endl;
        cout <<  "*******************************" << endl << endl;
        
        cout << "Pedigree: " << currentPed->filename() << endl;
        cout << "Phenotype: " << phenotype_filename << endl;
        cout << "Number of Subjects: " << eigenvalues.rows() << endl; 
        cout << "Total iterations: " << main_iteration_count << endl;
        cout << "Numerical Differentiation Delta: " << main_best_h << endl << endl;
      //  cout << "Final parameters\n";
        string parameter_names[6 + covariate_matrix.cols()*2];

        parameter_names[0] = trait_list[0] + "-h2r";
        int largest_name = parameter_names[0].length();

        parameter_names[1] = trait_list[1] + "-h2r";
        largest_name = ( parameter_names[1].length() > largest_name )  ? parameter_names[1].length() : largest_name;

        parameter_names[2] = trait_list[0] + "-sd";
        largest_name = ( parameter_names[2].length() > largest_name )  ? parameter_names[2].length() : largest_name;

        parameter_names[3] = trait_list[1] + "-sd";
        largest_name = ( parameter_names[3].length() > largest_name )  ? parameter_names[3].length() : largest_name;

        parameter_names[4] = "rhog";
        largest_name = ( parameter_names[4].length() > largest_name )  ? parameter_names[4].length() : largest_name;  

        parameter_names[5] = "rhoe";
        largest_name = ( parameter_names[5].length() > largest_name )  ? parameter_names[5].length() : largest_name;

        for(int i = 0 ; i < covariate_matrix.cols() - 1; i++){
            parameter_names[i + 6] = trait_list[0] + "-b" + cov_list[i];
            largest_name = ( parameter_names[6 + i].length() > largest_name )  ? parameter_names[6 + i].length() : largest_name;

            parameter_names[i + 6 + covariate_matrix.cols()] = trait_list[1] + "-b" + cov_list[i];
            largest_name = ( parameter_names[i + 6 + covariate_matrix.cols()].length() > largest_name )  ? parameter_names[i + 6 + covariate_matrix.cols()].length() : largest_name;            
        }

        parameter_names[6 + covariate_matrix.cols() - 1] = trait_list[0] + "-mean";
        largest_name = ( parameter_names[6 + covariate_matrix.cols() - 1].length() > largest_name )  ? parameter_names[6 + covariate_matrix.cols() - 1].length() : largest_name;


        parameter_names[6 + covariate_matrix.cols()*2 - 1] = trait_list[1] + "-mean";
        largest_name = ( parameter_names[6 + 2*covariate_matrix.cols() - 1].length() > largest_name )  ? parameter_names[6 + 2*covariate_matrix.cols() - 1].length() : largest_name;        
        largest_name = (14 > largest_name ) ? 14 : largest_name;
        cout << setw(largest_name + 4) << "Parameter"  << setw(20) << "Value" << setw(20) << "Standard Error" << endl << endl;
        
        cout << setw(largest_name + 4) << parameter_names[6 +covariate_matrix.cols() - 1] << setw(20) << final_beta(covariate_matrix.cols() - 1) << setw(20) << standard_errors(6 + covariate_matrix.cols() - 1) << endl;
        cout << setw(largest_name + 4) << parameter_names[6 + covariate_matrix.cols()*2 - 1] << setw(20) << final_beta(covariate_matrix.cols()*2  - 1)  << setw(20) << standard_errors(6 + 2*covariate_matrix.cols() - 1) << endl;
        for(int i = 0 ; i < covariate_matrix.cols() - 1; i++){
            cout << setw(largest_name + 4) << parameter_names[6 + i] << setw(20) << final_beta(i) << setw(20) << standard_errors(6 + i) << endl;
            cout << setw(largest_name + 4) << parameter_names[6 + covariate_matrix.cols() + i] << setw(20) << final_beta(i + covariate_matrix.cols()) << setw(20) << standard_errors(6 + covariate_matrix.cols() + i) << endl;
        }

        for(int i = 0 ; i < 6; i++){
            cout << setw(largest_name + 4) << parameter_names[i] << setw(20) << final_parameters(2 + i)<< setw(20) << standard_errors(i) << endl;
             
        }
        cout << setw(largest_name + 4) << "rhop" << setw(20) << rhop << endl << endl;
        cout << setw(largest_name + 4) << "loglik" << setw(20) << main_loglik << endl;
        cout << setw(largest_name + 4) << "rhog loglik" << setw(20) << rhog_loglik << endl;
        cout << setw(largest_name + 4) << "rhog p-value" << setw(20) << rhog_pvalue << endl;
        cout << setw(largest_name + 4) << "rhoe loglik" << setw(20) << rhoe_loglik << endl;
        cout << setw(largest_name + 4) << "rhoe p-value" << setw(20) << rhoe_pvalue << endl;
        
        


        string output_filename = "gen_corr-" + trait_list[0] + "-" + trait_list[1] + ".out";
        ofstream output_stream(output_filename);
        output_stream << "Parameter,Value,Standard Error\n";
        output_stream << parameter_names[6 +covariate_matrix.cols() - 1] << "," << final_beta(covariate_matrix.cols() - 1) << "," << standard_errors(6 + covariate_matrix.cols() - 1) << endl;
        output_stream << parameter_names[6 + covariate_matrix.cols()*2  - 1] << "," << final_beta(covariate_matrix.cols()*2  - 1) << "," << standard_errors(6 + 2*covariate_matrix.cols() - 1) << endl;
        for(int i = 0 ; i < covariate_matrix.cols() - 1; i++){
            output_stream << parameter_names[6 + i] << "," << final_beta(i) << ","<< standard_errors(6 + i) << endl;
            output_stream << parameter_names[6 + covariate_matrix.cols() + i] << "," << final_beta(i + covariate_matrix.cols()) << ","<< standard_errors(6 + covariate_matrix.cols() + i) << endl;
        }

        for(int i = 0 ; i < 6; i++){
            output_stream << parameter_names[i] << "," << final_parameters(i)<< ","<< standard_errors(i) << endl;
             
        } 
        output_stream << "rhop," << rhop << ",\n";
        output_stream << "loglik," << main_loglik << ",\n";
        output_stream << "rhog p-value," << rhog_pvalue << ",\n";
        output_stream << "rhog loglik," << rhog_loglik << ",\n";
        output_stream << "rhoe p-value," << rhoe_pvalue << ",\n"; 
        output_stream << "rhoe loglik," << rhoe_loglik << ",\n"; 
        output_stream.close();    
        
    }else{
        const double rhop = calculate_rho(final_parameters(4))*sqrt(calculate_constraint(final_parameters(0)))*\
             sqrt(calculate_constraint(final_parameters(1))) + calculate_rho(final_parameters(5))*\
             sqrt(1.0-calculate_constraint(final_parameters(0)))*sqrt(1.0-calculate_constraint(final_parameters(1)));
        final_parameters(0) = calculate_constraint(final_parameters(0));
        final_parameters(1) = calculate_constraint(final_parameters(1)); 
        final_parameters(4) = calculate_rho(final_parameters(4));
        final_parameters(5) = calculate_rho(final_parameters(5));  
        string rhop_command = "global gen_corr_RHOP ; set gen_corr_RHOP " + to_string(rhop);
        Solar_Eval(interp, rhop_command.c_str());

        cout << endl;
        cout << "*******************************" << endl;
        cout <<  "*     Genetic Correlation     *" << endl;
        cout <<  "*******************************" << endl << endl;
        
        cout << "Pedigree: " << currentPed->filename() << endl;
        cout << "Phenotype: " << phenotype_filename << endl;
        cout << "Number of Subjects: " << eigenvalues.rows() << endl; 
        cout << "Total iterations: " << main_iteration_count << endl;
        cout << "Numerical Differentiation Delta: " << main_best_h << endl << endl;
        
        string parameter_names[6 + covariate_matrix.cols()*2];

        parameter_names[0] = trait_list[0] + "-h2r";
        int largest_name = parameter_names[0].length();

        parameter_names[1] = trait_list[1] + "-h2r";
        largest_name = ( parameter_names[1].length() > largest_name )  ? parameter_names[1].length() : largest_name;

        parameter_names[2] = trait_list[0] + "-sd";
        largest_name = ( parameter_names[2].length() > largest_name )  ? parameter_names[2].length() : largest_name;

        parameter_names[3] = trait_list[1] + "-sd";
        largest_name = ( parameter_names[3].length() > largest_name )  ? parameter_names[3].length() : largest_name;

        parameter_names[4] = "rhog";
        largest_name = ( parameter_names[4].length() > largest_name )  ? parameter_names[4].length() : largest_name;  

        parameter_names[5] = "rhoe";
        largest_name = ( parameter_names[5].length() > largest_name )  ? parameter_names[5].length() : largest_name;

        for(int i = 0 ; i < covariate_matrix.cols() - 1; i++){
            parameter_names[i + 6] = trait_list[0] + "-b" + cov_list[i];
            largest_name = ( parameter_names[6 + i].length() > largest_name )  ? parameter_names[6 + i].length() : largest_name;

            parameter_names[i + 6 + covariate_matrix.cols()] = trait_list[1] + "-b" + cov_list[i];
            largest_name = ( parameter_names[i + 6 + covariate_matrix.cols()].length() > largest_name )  ? parameter_names[i + 6 + covariate_matrix.cols()].length() : largest_name;            
        }

        parameter_names[6 + covariate_matrix.cols() - 1] = trait_list[0] + "-mean";
        largest_name = ( parameter_names[6 + covariate_matrix.cols() - 1].length() > largest_name )  ? parameter_names[6 + covariate_matrix.cols() - 1].length() : largest_name;


        parameter_names[6 + covariate_matrix.cols()*2 - 1] = trait_list[1] + "-mean";
        largest_name = ( parameter_names[6 + 2*covariate_matrix.cols() - 1].length() > largest_name )  ? parameter_names[6 + 2*covariate_matrix.cols() - 1].length() : largest_name;        
        largest_name = (14 > largest_name ) ? 14 : largest_name;
       // cout << setw((largest_name + 4)*2) << "Final Parameters" << endl << endl;
        cout << setw(largest_name + 4) << "Parameter"  << setw(20) << "Value" << setw(20) << "Standard Error" << endl << endl;
        
        cout << setw(largest_name + 4) << parameter_names[6 +covariate_matrix.cols() - 1] << setw(20) << final_beta(covariate_matrix.cols() - 1) << setw(20) << standard_errors(6 + covariate_matrix.cols() - 1) << endl;
        cout << setw(largest_name + 4) << parameter_names[6 + 2*covariate_matrix.cols()  - 1] << setw(20) << final_beta(covariate_matrix.cols()*2  - 1) << setw(20) << standard_errors(6 + 2*covariate_matrix.cols() - 1) << endl;
        for(int i = 0 ; i < covariate_matrix.cols() - 1; i++){
            cout << setw(largest_name + 4) << parameter_names[6 + i] << setw(20) << final_beta(i) << setw(20) << standard_errors(6 + i) << endl;
            cout << setw(largest_name + 4) << parameter_names[6 + covariate_matrix.cols() + i] << setw(20) << final_beta(i + covariate_matrix.cols()) << setw(20) << standard_errors(6 + covariate_matrix.cols() + i) << endl;
        }

        for(int i = 0 ; i < 6; i++){
            cout << setw(largest_name + 4) << parameter_names[i] << setw(20) << final_parameters(i)<< setw(20) << standard_errors(i) << endl;
             
        }
        cout << setw(largest_name + 4) << "rhop" << setw(20) << rhop << endl;
        cout << setw(largest_name + 4) << "loglik" << setw(20) << main_loglik << endl;
 


        string output_filename = "gen_corr-" + trait_list[0] + "-" + trait_list[1] + ".out";
        ofstream output_stream(output_filename);
        output_stream << "Parameter,Value,Standard Error\n";
        output_stream << parameter_names[6 +covariate_matrix.cols() - 1] << "," << final_beta(covariate_matrix.cols() - 1)  << "," << standard_errors(6 + covariate_matrix.cols() - 1) << endl;
        output_stream << parameter_names[6 + covariate_matrix.cols()*2 - 1] << "," << final_beta(covariate_matrix.cols()*2  - 1)  << "," << standard_errors(6 + 2*covariate_matrix.cols() - 1) << endl;
        for(int i = 0 ; i < covariate_matrix.cols() - 1; i++){
            output_stream << parameter_names[6 + i] << "," << final_beta(i) << "," << standard_errors(6 + i) << endl;
            output_stream << parameter_names[6 + covariate_matrix.cols() + i] << "," << final_beta(i + covariate_matrix.cols()) << ","<< standard_errors(6 + covariate_matrix.cols() + i) << endl;
        }

        for(int i = 0 ; i < 6; i++){
            output_stream << parameter_names[i] << "," << final_parameters(i)<< ","<< standard_errors(i) << endl;
             
        } 
        output_stream << "rhop," << rhop << ",\n";
        output_stream << "loglik," << main_loglik << ",\n";
        output_stream.close();             

        }                        
                                
       
  
    delete file_data; 
    const char *error = 0;
    return error;               
} 
static void print_genetic_correlation_help(Tcl_Interp * interp){
    Solar_Eval(interp, "help gen_corr");
} 

extern "C" int run_genetic_correlation(ClientData clientData, Tcl_Interp * interp,
                          int argc, const char * argv[]){
  
    bool debug_mode = false;
 


    const char * evd_data_filename = 0;
    

    double h = 0.01;
    bool get_pvalues = false;
    for(int arg = 1 ;arg < argc ; arg++){
        if(!StringCmp(argv[arg], "help", case_ins) || !StringCmp(argv[arg], "-help", case_ins) || !StringCmp(argv[arg], "--help", case_ins)
           || !StringCmp(argv[arg], "h", case_ins) || !StringCmp(argv[arg], "-h", case_ins) || !StringCmp(argv[arg], "--help", case_ins)){
            print_genetic_correlation_help(interp);
            return TCL_OK;
        }else if (!StringCmp(argv[arg], "-d", case_ins) || !StringCmp(argv[arg], "--d", case_ins) || !StringCmp(argv[arg], "-debug", case_ins) ||
                  !StringCmp(argv[arg], "--debug", case_ins)){
            
            debug_mode  = true;
        }else if (!StringCmp(argv[arg], "-pvalues", case_ins) ||
                  !StringCmp(argv[arg], "--pvalues", case_ins)){
            
            get_pvalues  = true;
        }else if ((!StringCmp(argv[arg], "-evd_data", case_ins) || !StringCmp(argv[arg], "--evd_data", case_ins)) && arg + 1 < argc){
            evd_data_filename = argv[++arg];
        }else if ((!StringCmp(argv[arg], "-delta", case_ins) || !StringCmp(argv[arg], "--delta", case_ins)) && arg + 1 < argc){
            //evd_data_filename = argv[++arg];
            h = atof(argv[++arg]);
        }else{
            RESULT_LIT("Invalid argument enter see help");
            return TCL_ERROR;
        }
    }
    
    if(!loadedPed()){
        RESULT_LIT("No pedigree has been loaded");
        return TCL_ERROR;
    }
    const char * phenotype_filename = 0;
    phenotype_filename = Phenotypes::filenames();
    if(!phenotype_filename){
	RESULT_LIT("No phenotype file has bee loaded");
	return TCL_ERROR;
    }
    const char * pedigree_filename = currentPed->filename();
    
    if (Trait::Number_Of() != 2){
        RESULT_LIT( "Genetic correlation command requires two traits");
        return TCL_ERROR;
    }
    vector<string> trait_list;
    trait_list.push_back(string(Trait::Name(0)));
    trait_list.push_back(string(Trait::Name(1)));

    const char * error_msg = calculate_genetic_correlation(interp, trait_list,  phenotype_filename, get_pvalues, debug_mode, h, evd_data_filename);
    if(error_msg){
        RESULT_BUF(error_msg);
        return TCL_ERROR;
    }
    return TCL_OK;
}
extern "C" int transfer_pedigree_filename(ClientData clientData, Tcl_Interp * interp,
                          int argc, const char * argv[]){
      if(!loadedPed()){
        RESULT_LIT(0);
        return TCL_ERROR;
    }
    const char * pedigree_filename = currentPed->filename();
    RESULT_BUF(pedigree_filename);
    return TCL_OK;
    
} 
