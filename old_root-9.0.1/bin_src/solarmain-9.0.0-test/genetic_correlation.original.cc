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
   // const double e_x = exp(-x);
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
static double calculate_loglikelihood_param_all(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two, Eigen::MatrixXd covariate_matrix, Eigen::VectorXd lambda, Eigen::VectorXd parameters){
                   
    Eigen::VectorXd omega_one_one;// = omega_1_1 = sd_one*sd_one*(h2r_one*lambda + (1.0-h2r_one)*Eigen::VectorXd::Ones(lambda.rows()));
    Eigen::VectorXd omega_one_two;
     Eigen::VectorXd omega_two_two;
     const double h2r_one = calculate_constraint(parameters(0));
     const double h2r_two = calculate_constraint(parameters(1));
     const double rhog =  calculate_rho(parameters(4));
     const  double rhoe =  calculate_rho(parameters(5));
     const double sd_one = fabs(parameters(2));
     const  double sd_two = fabs(parameters(3));
      omega_one_one =  sd_one*sd_one*(h2r_one*lambda + (1.0-h2r_one)*Eigen::VectorXd::Ones(lambda.rows()));
      omega_one_two = sd_one*sd_two*(lambda*sqrt(h2r_one)*sqrt(h2r_two)*rhog + rhoe*sqrt(1.0-h2r_one)*sqrt(1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
      omega_two_two = sd_two*sd_two*(h2r_two*lambda + (1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
    
   // Eigen::VectorXd omega_one_two = omega_1_2 = sd_one*parameters(5)*(lambda*sqrt(h2r_one)*sqrt(h2r_two)*abs(rhog) + absrhoe*sqrt(1.0-h2r_one)*sqrt(1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
   /// Eigen::VectorXd omega_two_two = omega_2_2 = parameters(5)*parameters(5)*(h2r_two*lambda + (1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
    Eigen::VectorXd beta_one(covariate_matrix.cols());
    Eigen::VectorXd beta_two(covariate_matrix.cols());
    for(int i = 0 ; i < covariate_matrix.cols(); i++){
        beta_one(i) = parameters(6 + i);
        beta_two(i) = parameters(6 + i + covariate_matrix.cols());
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

        

    return -0.5*(part_one + part_two);
}
static double calculate_grad(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two,\
                        Eigen::MatrixXd covariate_matrix, Eigen::VectorXd lambda, \
                        Eigen::VectorXd all_parameters, double h,const int index, const bool is_hessian = false){

    double delta;
    if(is_hessian){
        delta = pow(DP, 0.66667)*max(abs(all_parameters(index)), 1.0);
    }else{
        delta = DP*max(abs(all_parameters(index)), 1.0);
    }
    Eigen::VectorXd positive_parameters = all_parameters;
    positive_parameters(index) += delta;
    Eigen::VectorXd negative_parameters =all_parameters;
    negative_parameters(index) -= delta;
    return (calculate_loglikelihood_param_all(Y_one,  Y_two,  covariate_matrix, \
     lambda,  positive_parameters) - calculate_loglikelihood_param_all(Y_one,  \
     Y_two,  covariate_matrix, lambda,  negative_parameters))/(2.0*delta);
}
static double calculate_hessian(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two,\
                        Eigen::MatrixXd covariate_matrix, Eigen::VectorXd lambda, \
                        Eigen::VectorXd all_parameters, double h,const int index_one,\
                        const int index_two){
    Eigen::VectorXd positive_parameters = all_parameters;
    double delta = pow(DP, 0.66667)*max(abs(all_parameters(index_two)), 1.0);
    positive_parameters(index_two) += delta;
    
    return ((calculate_grad(Y_one, Y_two, covariate_matrix, \
        lambda, positive_parameters, h, index_one, true)) - calculate_grad(Y_one, Y_two, covariate_matrix, \
        lambda, all_parameters, h, index_one, true))/(delta);
}
static double calculate_loglikelihood_e2(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two, Eigen::VectorXd Ones, Eigen::VectorXd lambda, const double beta_one, const double beta_two, \
                        const double sd_one, const double sd_two, const double e2_one, const double e2_two, const double rhog, const double rhoe, double h2r_one = -1, double h2r_two = -1){
   
                  
    Eigen::VectorXd omega_one_one;// = omega_1_1 = sd_one*sd_one*(h2r_one*lambda + (1.0-h2r_one)*Eigen::VectorXd::Ones(lambda.rows()));
    Eigen::VectorXd omega_one_two;
     Eigen::VectorXd omega_two_two;
       if (h2r_one == -1) 
            h2r_one = 1.0 - e2_one;
       if (h2r_two == -1)
           h2r_two = 1.0 - e2_two;
         omega_one_one  = sd_one*sd_one*(h2r_one*lambda + (e2_one)*Eigen::VectorXd::Ones(lambda.rows()));
        omega_one_two  = sd_one*sd_two*(lambda*sqrt(h2r_one)*sqrt(h2r_two)*rhog + rhoe*sqrt(e2_one)*sqrt(e2_two)*Eigen::VectorXd::Ones(lambda.rows()));
        omega_two_two  = sd_two*sd_two*((h2r_two)*lambda + (e2_two)*Eigen::VectorXd::Ones(lambda.rows()));
   
   // Eigen::VectorXd omega_one_two = omega_1_2 = sd_one*sd_two*(lambda*sqrt(h2r_one)*sqrt(h2r_two)*abs(rhog) + absrhoe*sqrt(1.0-h2r_one)*sqrt(1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
   /// Eigen::VectorXd omega_two_two = omega_2_2 = sd_two*sd_two*(h2r_two*lambda + (1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));

    Eigen::VectorXd residual_one = Y_one - Ones*beta_one;
    Eigen::VectorXd residual_two = Y_two - Ones*beta_two;

     
    
    double part_one = 0.0;

   Eigen::VectorXd omega_det = omega_one_one.cwiseProduct(omega_two_two) - omega_one_two.cwiseAbs2();
    for(int i = 0; i < lambda.rows() ; i++){
        part_one += log(abs(omega_det(i)));
    }

    compute_blockwise_diagonal_inversions(omega_one_one, omega_one_two, omega_two_two);
   

    const double part_two = residual_one.dot(residual_one.cwiseProduct(omega_one_one) + residual_two.cwiseProduct(omega_one_two)) + residual_two.dot(residual_two.cwiseProduct(omega_two_two) + residual_one.cwiseProduct(omega_one_two));// + 2.0*residual_one.dot(omega_one_two.cwiseProduct(residual_two));
  //  if (show_details){
  //      cout << "part one: " << part_one <<  " part two: " << part_two <<  " loglik: " << -0.5*(part_one + part_two) << endl;
 //   }
        

    return -0.5*(part_one + part_two);
}

static double calculate_loglikelihood_param_two(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two, Eigen::MatrixXd covariate_matrix, Eigen::VectorXd lambda, Eigen::VectorXd parameters, Eigen::VectorXd beta, const bool use_constraints = true){
                   
    Eigen::VectorXd omega_one_one;// = omega_1_1 = sd_one*sd_one*(h2r_one*lambda + (1.0-h2r_one)*Eigen::VectorXd::Ones(lambda.rows()));
    Eigen::VectorXd omega_one_two;
     Eigen::VectorXd omega_two_two;
     const double h2r_one = (use_constraints) ? calculate_constraint(parameters(2)) : parameters(2);
     const double h2r_two = (use_constraints) ? calculate_constraint(parameters(3)) : parameters(3);
     const double rhog = (use_constraints) ? calculate_rho(parameters(6)) : parameters(6);
     const  double rhoe =  (use_constraints) ? calculate_rho(parameters(7)) : parameters(7);
     const double sd_one = fabs(parameters(4));
     const  double sd_two = fabs(parameters(5));
      omega_one_one =  sd_one*sd_one*(h2r_one*lambda + (1.0-h2r_one)*Eigen::VectorXd::Ones(lambda.rows()));
      omega_one_two = sd_one*sd_two*(lambda*sqrt(h2r_one)*sqrt(h2r_two)*rhog + rhoe*sqrt(1.0-h2r_one)*sqrt(1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
      omega_two_two = sd_two*sd_two*(h2r_two*lambda + (1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
    
   // Eigen::VectorXd omega_one_two = omega_1_2 = sd_one*parameters(5)*(lambda*sqrt(h2r_one)*sqrt(h2r_two)*abs(rhog) + absrhoe*sqrt(1.0-h2r_one)*sqrt(1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
   /// Eigen::VectorXd omega_two_two = omega_2_2 = parameters(5)*parameters(5)*(h2r_two*lambda + (1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
    Eigen::VectorXd beta_one(covariate_matrix.cols());
    Eigen::VectorXd beta_two(covariate_matrix.cols());
    for(int i = 0 ; i < covariate_matrix.cols(); i++){
        beta_one(i) = beta(i);
        beta_two(i) = beta(i + covariate_matrix.cols());
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
static double calculate_loglikelihood_param(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two, Eigen::VectorXd Ones, Eigen::VectorXd lambda, Eigen::VectorXd parameters){
                   
    Eigen::VectorXd omega_one_one;// = omega_1_1 = sd_one*sd_one*(h2r_one*lambda + (1.0-h2r_one)*Eigen::VectorXd::Ones(lambda.rows()));
    Eigen::VectorXd omega_one_two;
     Eigen::VectorXd omega_two_two;
     double h2r_one = calculate_constraint(parameters(2));
     double h2r_two = calculate_constraint(parameters(3));
     double rhog = calculate_rho(parameters(6));
      double rhoe = calculate_rho(parameters(7));

      omega_one_one =  parameters(4)*parameters(4)*(h2r_one*lambda + (1.0-h2r_one)*Eigen::VectorXd::Ones(lambda.rows()));
      omega_one_two = parameters(4)*parameters(5)*(lambda*sqrt(h2r_one)*sqrt(h2r_two)*rhog + rhoe*sqrt(1.0-h2r_one)*sqrt(1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
      omega_two_two = parameters(5)*parameters(5)*(h2r_two*lambda + (1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
    
   // Eigen::VectorXd omega_one_two = omega_1_2 = sd_one*parameters(5)*(lambda*sqrt(h2r_one)*sqrt(h2r_two)*abs(rhog) + absrhoe*sqrt(1.0-h2r_one)*sqrt(1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
   /// Eigen::VectorXd omega_two_two = omega_2_2 = parameters(5)*parameters(5)*(h2r_two*lambda + (1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));

    Eigen::VectorXd residual_one = Y_one - Ones*parameters(0);
    Eigen::VectorXd residual_two = Y_two - Ones*parameters(1);

     
    
    double part_one = 0.0;

   Eigen::VectorXd omega_det = omega_one_one.cwiseProduct(omega_two_two) - omega_one_two.cwiseAbs2();
    for(int i = 0; i < lambda.rows() ; i++){
        part_one += log(abs(omega_det(i)));
    }

    compute_blockwise_diagonal_inversions(omega_one_one, omega_one_two, omega_two_two);
   

    const double part_two = residual_one.dot(residual_one.cwiseProduct(omega_one_one)) + residual_two.dot(residual_two.cwiseProduct(omega_two_two)) + 2.0*residual_one.dot(omega_one_two.cwiseProduct(residual_two));

        

    return -0.5*(part_one + part_two);
}
static double calculate_gradient_param(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two, Eigen::VectorXd Ones, Eigen::VectorXd lambda, Eigen::VectorXd parameters, int index, double h){
     Eigen::VectorXd positive_parameters = parameters;
     positive_parameters(index) += h;
      Eigen::VectorXd negative_parameters = parameters;
      negative_parameters(index) -= h;
    return (calculate_loglikelihood_param(Y_one, Y_two, Ones, lambda, positive_parameters)  - calculate_loglikelihood_param(Y_one, Y_two, Ones, lambda, negative_parameters))/(2.0*h) ;  
}  
static double calculate_loglikelihood(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two, Eigen::VectorXd Ones, Eigen::VectorXd lambda, const double beta_one, const double beta_two, \
                        const double sd_one, const double sd_two, const double h2r_one, const double h2r_two, const double rhog, const double rhoe, bool is_hessian = false, bool show_details = false){
                  
    Eigen::VectorXd omega_one_one;// = omega_1_1 = sd_one*sd_one*(h2r_one*lambda + (1.0-h2r_one)*Eigen::VectorXd::Ones(lambda.rows()));
    Eigen::VectorXd omega_one_two;
     Eigen::VectorXd omega_two_two;

     omega_one_one = sd_one*sd_one*(h2r_one*lambda + (1.0-h2r_one)*Eigen::VectorXd::Ones(lambda.rows()));
     omega_one_two  = sd_one*sd_two*(lambda*sqrt(h2r_one)*sqrt(h2r_two)*rhog + rhoe*sqrt(1.0-h2r_one)*sqrt(1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
     omega_two_two  = sd_two*sd_two*(h2r_two*lambda + (1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
    
   // Eigen::VectorXd omega_one_two = omega_1_2 = sd_one*sd_two*(lambda*sqrt(h2r_one)*sqrt(h2r_two)*abs(rhog) + absrhoe*sqrt(1.0-h2r_one)*sqrt(1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));
   /// Eigen::VectorXd omega_two_two = omega_2_2 = sd_two*sd_two*(h2r_two*lambda + (1.0-h2r_two)*Eigen::VectorXd::Ones(lambda.rows()));

    Eigen::VectorXd residual_one = Y_one - Ones*beta_one;
    Eigen::VectorXd residual_two = Y_two - Ones*beta_two;

     
    
    double part_one = 0.0;

   Eigen::VectorXd omega_det = omega_one_one.cwiseProduct(omega_two_two) - omega_one_two.cwiseAbs2();
    for(int i = 0; i < lambda.rows() ; i++){
        part_one += log(abs(omega_det(i)));
    }

    compute_blockwise_diagonal_inversions(omega_one_one, omega_one_two, omega_two_two);
   

    const double part_two = residual_one.dot(residual_one.cwiseProduct(omega_one_one) + residual_two.cwiseProduct(omega_one_two)) + residual_two.dot(residual_two.cwiseProduct(omega_two_two) + residual_one.cwiseProduct(omega_one_two));// + 2.0*residual_one.dot(omega_one_two.cwiseProduct(residual_two));
    if (show_details){
        cout << "part one: " << part_one <<  " part two: " << part_two <<  " loglik: " << -0.5*(part_one + part_two) << endl;
    }
        

    return -0.5*(part_one + part_two);
}
static Eigen::VectorXd calculate_errors(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two, Eigen::VectorXd Ones, Eigen::VectorXd lambda, double mu_one, double mu_two,\
                                        double t_one, double t_two, double sd_one, double sd_two, double i_one, double i_two, double h){
    double mu_one_h = h;
    double mu_two_h = h;
    double sd_one_h = h;
    double sd_two_h = h;
    double t_one_h = h;
    double t_two_h = h;
    double i_one_h = h;
    double i_two_h = h;
    double e2_one = 1.0 - calculate_constraint(t_one);
    double e2_two = 1.0 - calculate_constraint(t_two);
    double e2_one_h = h;
    double e2_two_h = h;
    double h2r_one_h = h;
    double h2r_two_h = h;
    double rhog_h = h;
    double rhoe_h = h;
    double h2r_one = 1.0 - e2_one;
    double h2r_two = 1.0 - e2_two; 
    double rhog = calculate_rho(i_one);
    double rhoe = calculate_rho(i_two);     
    if(e2_one + h > 1.0){
        e2_one_h = 1.0-e2_one;
    }
    if(e2_one - h < 0.0){
        e2_one_h = e2_one;
    }
    
    if(e2_two + h > 1.0){
        e2_two_h = 1.0-e2_two;
    }
    if(e2_two - h < 0.0){
        e2_two_h = e2_two;
    }   
    
    if(h2r_one + h > 1.0){
        h2r_one_h = 1.0-h2r_one;
    }
    if(h2r_one - h < 0.0){
        h2r_one_h = h2r_one;
    }
    
    if(h2r_two + h > 1.0){
        h2r_two_h = 1.0-h2r_two;
    }
    if(h2r_two - h < 0.0){
        h2r_two_h = h2r_two;
    }      

    if(rhog + h > 1.0){
        rhog_h = 1.0-rhog;
    }
    if(rhog - h < -1.0){
        rhog_h = rhog + 1.0;
    } 
     
    if(rhoe + h > 1.0){
        rhoe_h = 1.0-rhoe;
    }
    if(rhoe - h < -1.0){
        rhoe_h = rhoe + 1.0;
    } 
    Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(10, 10);
    double loglik  = calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, rhog, rhoe);

     hessian(0, 0) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)))/(mu_one_h*mu_one_h);
    
        hessian(1, 1) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two + h, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two - h, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)))/(mu_two_h*mu_two_h);

    hessian(2, 2) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two, sd_one, sd_two, e2_one + e2_one_h, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one - e2_one_h, e2_two, calculate_rho(i_one), calculate_rho(i_two)))/(e2_one_h*e2_one_h);
    
    hessian(3, 3) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two, sd_one, sd_two, e2_one, e2_two + e2_two_h, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two - e2_two_h, calculate_rho(i_one), calculate_rho(i_two)))/(e2_two_h*e2_two_h);             

    hessian(4, 4) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two, sd_one + sd_one_h, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one - sd_one_h, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)))/(sd_one_h*sd_one_h);

    hessian(5, 5) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two, sd_one, sd_two + sd_two_h, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two - sd_two_h, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)))/(sd_two_h*sd_two_h);

    hessian(6, 6) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two, sd_one, sd_two, e2_one, e2_two, rhog + rhog_h, calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, rhog - rhog_h, calculate_rho(i_two)))/(rhog_h*rhog_h); 

    hessian(7, 7) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), rhoe + rhoe_h) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), rhoe - rhoe_h))/(rhoe_h*rhoe_h);
    
    hessian(0, 1) = hessian(1, 0) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two + mu_two_h, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two - mu_two_h, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)))/((2.0*mu_one_h*mu_two_h) + 0.5*(hessian(0,0) + hessian(1, 1)));         
        

    hessian(0, 2) = hessian(2, 0) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two , sd_one, sd_two, e2_one + e2_one_h, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik  + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two, e2_one - e2_one_h, calculate_constraint(t_two ), calculate_rho(i_one), calculate_rho(i_two)))/((2.0*e2_one_h*mu_one_h) + 0.5*(hessian(0, 0) + hessian(2, 2)));


    hessian(0, 3) = hessian(3, 0) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two , sd_one, sd_two, e2_one, e2_two + e2_two_h, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two, e2_one, e2_two - e2_two_h, calculate_rho(i_one), calculate_rho(i_two)))/((2.0*e2_two_h*mu_one_h) + 0.5*(hessian(0,0) + hessian(3, 3)));    
    
 

    hessian(0, 4) = hessian(4, 0) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two , sd_one + sd_one_h, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one -sd_one_h, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)))/((2.0*sd_one_h*mu_one_h) + 0.5*(hessian(0, 0) + hessian(4, 4)));
    

  
  
      hessian(0, 5) = hessian(5, 0) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two , sd_one, sd_two + sd_two_h, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two -sd_two_h, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)))/((2.0*sd_two_h*mu_one_h) + 0.5*(hessian(0, 0) + hessian(5, 5))); 
 

 
     hessian(0, 6) = hessian(6, 0) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two , sd_one, sd_two, e2_one, e2_two, rhog  + rhog_h, calculate_rho(i_two)) - 2.0*loglik  + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two, e2_one, e2_two, rhog -  rhog_h, calculate_rho(i_two)))/((2.0*rhog_h*mu_one_h) + 0.5*(hessian(0, 0) + hessian(6, 6)));  


     hessian(0, 7) = hessian(7, 0) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two , sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), rhoe + rhoe_h) - 2.0*loglik  + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one ), rhoe - rhoe_h))/((2.0*rhoe_h*mu_one_h) + 0.5*(hessian(0,0) + hessian(7, 7)));     
    


    hessian(1, 2) = hessian(2, 1) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two + mu_two_h, sd_one, sd_two, e2_one + e2_one_h, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik  + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two - mu_two_h, sd_one, sd_two, e2_one - e2_one_h, e2_two, calculate_rho(i_one), calculate_rho(i_two)))/((2.0*e2_one_h*mu_two_h) + 0.5*(hessian(1, 1) + hessian(2, 2)));       
    


    hessian(1, 3) = hessian(3, 1) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two + mu_two_h, sd_one, sd_two, e2_one, e2_two + e2_two_h, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik  + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two - mu_two_h, sd_one, sd_two, e2_one, e2_two - e2_two_h, calculate_rho(i_one), calculate_rho(i_two)))/((2.0*e2_two_h*mu_two_h) + 0.5*(hessian(1, 1) + hessian(3, 3)));   
    

 
    hessian(1, 4) = hessian(4, 1) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two + mu_two_h, sd_one + sd_one_h, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two - mu_two_h, sd_one -sd_one_h, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)))/((2.0*sd_one_h*mu_two_h) + 0.5*(hessian(1, 1) + hessian(4, 4))); 
    

    hessian(1, 5) = hessian(5, 1) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two + mu_two_h , sd_one, sd_two + sd_two_h, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik  + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two - mu_two_h, sd_one, sd_two -sd_two_h, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)))/((2.0*sd_two_h*mu_two_h) + (hessian(1, 1) + hessian(5, 5)));   
 

 
 
     hessian(1, 6) = hessian(6, 1) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two + mu_two_h, sd_one, sd_two, e2_one, e2_two, rhog  + rhog_h, calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two - mu_two_h, sd_one, sd_two, e2_one, e2_two, rhog -  rhog_h, calculate_rho(i_two)))/((2.0*rhog_h*mu_two_h) + 0.5*(hessian(1, 1) + hessian(6, 6)));    


     hessian(1, 7) = hessian(7, 1) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two + mu_two_h , sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), rhoe + rhoe_h) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two - mu_two_h, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one ), rhoe - rhoe_h))/((2.0*rhoe_h*mu_two_h) + 0.5*(hessian(1, 1) + hessian(7, 7)));   
    
     
    
    hessian(2, 3) = hessian(3, 2) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two, e2_one + e2_one_h, e2_two + e2_two_h, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one - e2_one_h, e2_two - e2_two_h, calculate_rho(i_one), calculate_rho(i_two)))/((2.0*e2_one_h*e2_two_h) + 0.5*(hessian(2, 2) + hessian(3, 3)));
    



    hessian(2, 4) = hessian(4, 2) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one + sd_one_h, sd_two, e2_one + e2_one_h, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik  + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one - sd_one_h, sd_two, e2_one - e2_one_h, e2_two, calculate_rho(i_one), calculate_rho(i_two)))/((2.0*e2_one_h*sd_one_h) + 0.5*(hessian(2, 2) + hessian(4, 4)));
    

  
      hessian(2, 5) = hessian(5, 2) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two + sd_two_h, e2_one + h, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two - sd_two_h, e2_one - e2_one_h, e2_two, calculate_rho(i_one), calculate_rho(i_two)))/((2.0*e2_one_h*sd_two_h) + 0.5*(hessian(2, 2) + hessian(5, 5)));    
    

  
      hessian(2, 6) = hessian(6, 2) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two, e2_one + e2_one_h, e2_two, rhog  + rhog_h, calculate_rho(i_two)) - 2.0*loglik  + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one - e2_one_h, e2_two, rhog -  rhog_h, calculate_rho(i_two)))/((2.0*e2_one_h*rhog_h) + 0.5*(hessian(2, 2) + hessian(6, 6)));
    

    hessian(2, 7) = hessian(7, 2) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two, e2_one + e2_one_h, e2_two, calculate_rho(i_one), rhoe + rhoe_h) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one - e2_one_h, e2_two, calculate_rho(i_one), rhoe - rhoe_h))/((2.0*e2_one_h*rhoe_h) + 0.5*(hessian(2, 2) + hessian(7, 7))); 
       
    

    hessian(3, 4) = hessian(4, 3) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one + sd_one_h, sd_two, e2_one, e2_two + e2_two_h, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one - sd_one_h, sd_two, e2_one, e2_two - e2_two_h, calculate_rho(i_one), calculate_rho(i_two)))/((2.0*e2_two_h*sd_one_h) + 0.5*(hessian(3, 3) + hessian(4, 4)));
    


    hessian(3, 5) = hessian(5, 3) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two + sd_two_h, e2_one, e2_two + e2_two_h, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two - sd_two_h, e2_one, e2_two - e2_two_h, calculate_rho(i_one), calculate_rho(i_two)))/((2.0*t_two_h*sd_two_h) + 0.5*(hessian(3, 3) + hessian(5, 5)));
    

     hessian(3, 6) = hessian(6, 3) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two, e2_one, e2_two + e2_two_h, rhog  + rhog_h, calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two - e2_two_h, rhog -  rhog_h, calculate_rho(i_two)))/((2.0*e2_two_h*rhog_h) + 0.5*(hessian(3, 3) + hessian(6, 6))); 
    


    hessian(3, 7) = hessian(7, 3) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two, e2_one, e2_two + e2_two_h, calculate_rho(i_one), rhoe + rhoe_h) - 2.0*loglik  + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two - e2_two_h, calculate_rho(i_one), rhoe - rhoe_h))/((2.0*e2_two_h*rhoe_h) + 0.5*(hessian(3, 3) + hessian(7, 7)));
    


    hessian(4, 5) = hessian(5, 4) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one + sd_one_h, sd_two + sd_two_h, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one - sd_one_h, sd_two - sd_two_h, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two)))/((2.0*sd_one_h*sd_two_h) + 0.5*(hessian(4, 4) + hessian(5, 5))); 
    

     hessian(4, 6) = hessian(6, 4) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one + sd_one_h, sd_two, e2_one, e2_two, rhog  + rhog_h, calculate_rho(i_two)) - 2.0*loglik  + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one - sd_one_h, sd_two, e2_one, e2_two, rhog - rhog_h, calculate_rho(i_two)))/((2.0*sd_one_h*rhog_h) + 0.5*(hessian(4, 4) + hessian(6, 6))); 
    


    hessian(4, 7) = hessian(7, 4) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one + sd_one_h, sd_two, e2_one, e2_two, calculate_rho(i_one), rhoe + rhoe_h) - 2.0*loglik  + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one - sd_one_h, sd_two, e2_one, e2_two, calculate_rho(i_one), rhoe - rhoe_h))/((2.0*sd_one_h*rhoe_h) + 0.5*(hessian(4, 4) + hessian(7,7))); 
    

     hessian(5, 6) = hessian(6, 5) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two + sd_two_h, e2_one, e2_two, rhog  + rhog_h, calculate_rho(i_two)) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two  - sd_two_h, e2_one, e2_two, rhog - rhog_h, calculate_rho(i_two)))/((2.0*sd_two_h*rhog_h) + 0.5*(hessian(5, 5) + hessian(6, 6)));    



     hessian(5, 7) = hessian(7, 5) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one , sd_two + sd_two_h, e2_one, e2_two, calculate_rho(i_one), rhoe + rhoe_h) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two - sd_two_h, e2_one, e2_two, calculate_rho(i_one), rhoe - rhoe_h))/((2.0*sd_two_h*rhoe_h) + 0.5*(hessian(5, 5) + hessian(7, 7)));  


    hessian(6, 7) = hessian(7, 6) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one , sd_two, e2_one, e2_two, rhog  + rhog_h, rhoe + rhoe_h) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, rhog - rhog_h, rhoe - rhoe_h))/((2.0*rhoe_h*rhog_h) + 0.5*(hessian(6,6) + hessian(7, 7)));
    
    hessian(8, 8) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one + h2r_one_h) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one - h2r_one_h))/(h2r_one_h*h2r_one_h);
    
    hessian(9, 9) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one, h2r_two + h2r_two_h) - 2.0*loglik + calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one, h2r_two - h2r_two_h))/(h2r_two_h*h2r_two_h);   

    hessian(8, 0) = hessian(0, 8) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one + h2r_one_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one - h2r_one_h))/((2.0*mu_one_h*h2r_one_h) + 0.5*(hessian(0,0) + hessian(8, 8)));

    hessian(8, 1) = hessian(1, 8) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two + h, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one + h2r_one_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two - h, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one - h2r_one_h))/((2.0*h2r_one_h*mu_two_h) + 0.5*(hessian(1,1) + hessian(8, 8)));
  
    hessian(8, 2) = hessian(2, 8) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one + e2_one_h, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one + h2r_one_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one - e2_one_h, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one - h2r_one_h))/((2.0*h2r_one_h*e2_one_h) + 0.5*(hessian(2,2) + hessian(8, 8)));    
    hessian(8, 3) = hessian(3, 8) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two + e2_two_h, calculate_rho(i_one), calculate_rho(i_two), h2r_one + h2r_one_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two - e2_two_h, calculate_rho(i_one), calculate_rho(i_two), h2r_one - h2r_one_h))/((2.0*h2r_one_h*e2_two_h) + 0.5*(hessian(3,3) + hessian(8, 8))); 
    
    hessian(8, 4) = hessian(4, 8) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one + h, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one + h2r_one_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one-h, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one - h2r_one_h))/((2.0*mu_one_h*h2r_one_h) + 0.5*(hessian(4,4) + hessian(8, 8)));  
    
    hessian(8, 5) = hessian(5, 8) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two + h, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one + h2r_one_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two - h, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one - h2r_one_h))/((2.0*mu_one_h*h2r_one_h) + 0.5*(hessian(5,5) + hessian(8, 8))); 

    hessian(8, 6) = hessian(6, 8) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, rhog  + rhog_h, calculate_rho(i_two), h2r_one + h2r_one_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, rhog - rhog_h, calculate_rho(i_two), h2r_one - h2r_one_h))/((2.0*h2r_one_h*rhog_h) + 0.5*(hessian(6,6) + hessian(8, 8)));  
    
    hessian(8, 7) = hessian(7, 8) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, rhog, rhoe + rhoe_h, h2r_one + h2r_one_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, rhog, rhoe-rhoe_h, h2r_one - h2r_one_h))/((2.0*rhoe_h*h2r_one_h) + 0.5*(hessian(7,7) + hessian(8, 8))); 
    
    hessian(8, 9) = hessian(9, 8) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, rhog, rhoe, h2r_one + h2r_one_h, h2r_two + h2r_two_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two , e2_one, e2_two, rhog, rhoe, h2r_one - h2r_one_h, h2r_two - h2r_two_h))/((2.0*h2r_one_h*h2r_two_h) + 0.5*(hessian(9,9) + hessian(8, 8))); 
    

    hessian(9, 0) = hessian(0, 9) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one, h2r_two + h2r_two_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one,  h2r_two - h2r_two_h))/((2.0*mu_one_h*h2r_two_h) + 0.5*(hessian(0,0) + hessian(9, 9)));

    hessian(9, 1) = hessian(1, 9) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two + h, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one,  h2r_two + h2r_two_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two - h, sd_one, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one,  h2r_two - h2r_two_h))/((2.0*mu_one_h*h2r_two_h) + 0.5*(hessian(1,1) + hessian(9, 9)));
  
    hessian(9, 2) = hessian(2, 9) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one + e2_one_h, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one,  h2r_two + h2r_two_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one - e2_one_h, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one,  h2r_two - h2r_two_h))/((2.0*e2_one_h*h2r_two_h) + 0.5*(hessian(2,2) + hessian(9, 9)));    
    hessian(9, 3) = hessian(3, 9) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two + e2_two_h, calculate_rho(i_one), calculate_rho(i_two), h2r_one,  h2r_two + h2r_two_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two - e2_two_h, calculate_rho(i_one), calculate_rho(i_two), h2r_one,  h2r_two - h2r_two_h))/((2.0*e2_two_h*h2r_two_h) + 0.5*(hessian(3,3) + hessian(9, 9))); 
    
    hessian(9, 4) = hessian(4, 9) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one + h, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one,  h2r_two + h2r_two_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one-h, sd_two, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one,  h2r_two - h2r_two_h))/((2.0*mu_one_h*h2r_two_h) + 0.5*(hessian(4,4) + hessian(9, 9)));  
    
    hessian(9, 5) = hessian(5, 9) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two + h, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one,  h2r_two + h2r_two_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two - h, e2_one, e2_two, calculate_rho(i_one), calculate_rho(i_two), h2r_one,  h2r_two - h2r_two_h))/((2.0*mu_one_h*h2r_two_h) + 0.5*(hessian(5,5) + hessian(9, 9))); 

    hessian(9, 6) = hessian(6, 9) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, rhog  + rhog_h, calculate_rho(i_two), h2r_one,  h2r_two + h2r_two_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, rhog - rhog_h, calculate_rho(i_two), h2r_one,  h2r_two - h2r_two_h))/((2.0*rhog_h*h2r_two_h) + 0.5*(hessian(6,6) + hessian(9, 9)));  
    
    hessian(9, 7) = hessian(7, 9) = (calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, rhog, rhoe + rhoe_h, h2r_one,  h2r_two + h2r_two_h) - 2.0*loglik+ calculate_loglikelihood_e2(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, e2_one, e2_two, rhog, rhoe-rhoe_h, h2r_one,  h2r_two - h2r_two_h))/((2.0*rhoe_h*h2r_two_h) + 0.5*(hessian(7,7) + hessian(9, 9)));    
    
                    
    return (-hessian).inverse().diagonal().cwiseAbs().cwiseSqrt();
                                        
}
static void calculate_gradient(Eigen::VectorXd & gradient, Eigen::VectorXd Y_one, Eigen::VectorXd Y_two,\
                                                                Eigen::VectorXd Ones,  Eigen::VectorXd lambda, double mu_one, double mu_two, double t_one, double t_two,\
                                                                double sd_one, double sd_two, double i_one, double i_two, double h){
    gradient.resize(8);                                                           
    double h2r_one = calculate_constraint(t_one);  
    double h2r_two = calculate_constraint(t_two);   
    double rhog = calculate_rho(i_one);  
    double rhoe = calculate_rho(i_two);                                                   
    double mu_one_h = h;
    double mu_two_h = h;
    double t_one_h = h;
    double t_two_h = h;
    double sd_one_h = h;
    double sd_two_h = h;
    double i_one_h = h;
    double i_two_h = h;

    gradient(0) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)))/(2.0*mu_one_h) ;
                                                                        
    gradient(1) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two + mu_two_h, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two - mu_two_h, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)))/(2.0*mu_two_h) ;
                                                                 
    gradient(2) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, h2r_one + h, calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two,h2r_one - h, calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)))/(2.0*t_one_h) ; 
                                                                    
    gradient(3) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two , sd_one, sd_two, calculate_constraint(t_one), h2r_two + h, calculate_rho(i_one), calculate_rho(i_two)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one), h2r_two - h, calculate_rho(i_one), calculate_rho(i_two)))/(2.0*t_two_h) ;  
                                                                    
    gradient(4) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two , sd_one + sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one - sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)))/(2.0*sd_one_h) ; 
     
    gradient(5) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two , sd_one, sd_two + sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two - sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)))/(2.0*sd_two_h) ;    
    
    gradient(6) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two , sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), rhog + h, calculate_rho(i_two)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), rhog - h, calculate_rho(i_two)))/(2.0*i_one_h) ;  
 
    gradient(7) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two , sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), rhoe + h) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), rhoe - h))/(2.0*i_two_h) ; 
    
}
static void genetic_correlation_calculate_hessian_and_gradient_multithread_four(Eigen::MatrixXd & hessian, Eigen::VectorXd & gradient, Eigen::VectorXd Y_one, Eigen::VectorXd Y_two,\
                                                                Eigen::MatrixXd covariate_matrix,  Eigen::VectorXd lambda, Eigen::VectorXd parameters, Eigen::VectorXd beta){
                                                                
double loglik =  calculate_loglikelihood_param_two(Y_one,  Y_two,  covariate_matrix,  lambda,  parameters ,beta);
int index_map[6] = {2, 3, 4, 5, 6, 7};
double diff[6];
for(int i = 0 ; i < 6; i++){
    const double parameter = fabs(parameters(index_map[i]));
    diff[i] = DP*((parameter != 0.0 ) ? parameter : 1.0 );
    //diff[i] = DP*max(fabs(parameters(index_map[i])), 1.0);
}
//#pragma omp parallel for
    for(int index = 0; index < 6 ; index++){
        Eigen::VectorXd positive_parameters = parameters;
        positive_parameters(index_map[index]) += diff[index];
        Eigen::VectorXd negative_parameters = parameters;
        negative_parameters(index_map[index]) -= diff[index];
        double positive_loglik = calculate_loglikelihood_param_two(Y_one,  Y_two,  covariate_matrix,  lambda,  positive_parameters,beta);
        double negative_loglik = calculate_loglikelihood_param_two(Y_one,  Y_two,  covariate_matrix,  lambda,  negative_parameters,beta);
        gradient(index) =  (positive_loglik - negative_loglik) /(2.0*diff[index]);
        hessian(index, index) =  (positive_loglik - 2.0*loglik + negative_loglik) /(diff[index]*diff[index]);
    }
  
//#pragma omp parallel for
    for(int i = 0 ; i < 6 ; i++){
       /* const long v = 2*8+1;
        int j = floor( ( v - sqrt( (double)(v*v - 8*index) ) ) / 2.0 ) ;
        int i = index - 8*j + j*(j-1)/2  + j;
        if (j == i)
            continue;*/
        for(int j = i  + 1; j < 6 ; j++){
            Eigen::VectorXd positive_parameters = parameters;
            positive_parameters(index_map[i]) += diff[i];
            positive_parameters(index_map[j]) += diff[j];
            Eigen::VectorXd negative_parameters = parameters; 
            negative_parameters(index_map[i]) -= diff[i];
            negative_parameters(index_map[j]) -= diff[j];              
            double positive_loglik = calculate_loglikelihood_param_two(Y_one,  Y_two,  covariate_matrix,  lambda,  positive_parameters,beta);
            double negative_loglik = calculate_loglikelihood_param_two(Y_one,  Y_two,  covariate_matrix,  lambda,  negative_parameters,beta);
            hessian(i, j) = hessian(j, i) = (positive_loglik - 2.0*loglik + negative_loglik)/((2.0*diff[i]*diff[j]) + 0.5*(hessian(i, i) + hessian(j, j)));
        }
       
    }                                                                  
                                                                
} 
static void genetic_correlation_calculate_hessian_and_gradient_multithread_three(Eigen::MatrixXd & hessian, Eigen::VectorXd & gradient, Eigen::VectorXd Y_one, Eigen::VectorXd Y_two,\
                                                                Eigen::MatrixXd covariate_matrix,  Eigen::VectorXd lambda, Eigen::VectorXd parameters, Eigen::VectorXd beta, double h){
                                                                
//double loglik =  calculate_loglikelihood_param_two(Y_one,  Y_two,  covariate_matrix,  lambda,  parameters ,beta);
int index_map[6] = {2, 3, 4, 5, 6, 7};
Eigen::VectorXd all_parameters(6 + covariate_matrix.cols()*2);
for(int i = 0 ; i < 6; i++){
    all_parameters(i) = parameters(index_map[i]);
}
for(int i = 0; i < covariate_matrix.cols(); i++){
    all_parameters(6 + i) = beta(i);
    all_parameters(6 + i + covariate_matrix.cols()) = beta(i + covariate_matrix.cols());
}
int index_map_two[4] = {0,  1, 4, 5};
//#pragma omp parallel for
    for(int index = 0; index < 4 ; index++){
        //Eigen::VectorXd positive_parameters = parameters;
        //positive_parameters(index_map_two[index]) += h;
        //Eigen::VectorXd negative_parameters = parameters;
       // negative_parameters(index_map_two[index]) -= h;
       // double positive_loglik = calculate_loglikelihood_param_two(Y_one,  Y_two,  covariate_matrix,  lambda,  positive_parameters,beta);
       // double negative_loglik = calculate_loglikelihood_param_two(Y_one,  Y_two,  covariate_matrix,  lambda,  negative_parameters,beta);
        gradient(index) = calculate_grad(Y_one, Y_two,\
                        covariate_matrix, lambda, \
                        all_parameters,  h, index_map_two[index]);
      //  gradient(index) =  (positive_loglik - negative_loglik) /(2.0*h);
        hessian(index, index) = calculate_hessian(Y_one, Y_two,\
                        covariate_matrix, lambda, \
                        all_parameters,  h, index_map_two[index],\
                         index_map_two[index]);
        //hessian(index, index) =  (positive_loglik - 2.0*loglik + negative_loglik) /(h*h);
    }
  
//#pragma omp parallel for
    for(int i = 0 ; i < 4 ; i++){
       /* const long v = 2*8+1;
        int j = floor( ( v - sqrt( (double)(v*v - 8*index) ) ) / 2.0 ) ;
        int i = index - 8*j + j*(j-1)/2  + j;
        if (j == i)
            continue;*/
        for(int j = i  + 1; j < 4 ; j++){
          /*  Eigen::VectorXd positive_parameters = parameters;
            positive_parameters(index_map_two[i]) += h;
            positive_parameters(index_map_two[j]) += h;
            Eigen::VectorXd negative_parameters = parameters; 
            negative_parameters(index_map_two[i]) -= h;
            negative_parameters(index_map_two[j]) -= h;              
            double positive_loglik = calculate_loglikelihood_param_two(Y_one,  Y_two,  covariate_matrix,  lambda,  positive_parameters,beta);
            double negative_loglik = calculate_loglikelihood_param_two(Y_one,  Y_two,  covariate_matrix,  lambda,  negative_parameters,beta);
            hessian(i, j) = hessian(j, i) = (positive_loglik - 2.0*loglik + negative_loglik)/((2.0*h*h) + 0.5*(hessian(i, i) + hessian(j, j)));*/
         hessian(i, j) = hessian(j, i) = calculate_hessian(Y_one, Y_two,\
                        covariate_matrix, lambda, \
                        all_parameters,  h, index_map_two[i],\
                         index_map_two[j]);
        }
       
    }                                                                  
                                                                
} 


static void genetic_correlation_calculate_hessian_and_gradient_multithread_two(Eigen::MatrixXd & hessian, Eigen::VectorXd & gradient, Eigen::VectorXd Y_one, Eigen::VectorXd Y_two,\
                                                                Eigen::VectorXd Ones,  Eigen::VectorXd lambda, Eigen::VectorXd parameters, double h){
                                                                
double loglik =  calculate_loglikelihood_param(Y_one,  Y_two,  Ones,  lambda,  parameters);
int index_map[4] = {2, 3, 6, 7};
//#pragma omp parallel for
    for(int index = 0; index < 4 ; index++){
        Eigen::VectorXd positive_parameters = parameters;
        positive_parameters(index_map[index]) += h;
        Eigen::VectorXd negative_parameters = parameters;
        negative_parameters(index_map[index]) -= h;
        double positive_loglik = calculate_loglikelihood_param(Y_one,  Y_two,  Ones,  lambda,  positive_parameters);
        double negative_loglik = calculate_loglikelihood_param(Y_one,  Y_two,  Ones,  lambda,  negative_parameters);
        gradient(index) =  (positive_loglik - negative_loglik) /(2.0*h);
        hessian(index, index) =  (positive_loglik - 2.0*loglik + negative_loglik) /(h*h);
    }
  
//#pragma omp parallel for
    for(int i = 0 ; i < 4 ; i++){
       /* const long v = 2*8+1;
        int j = floor( ( v - sqrt( (double)(v*v - 8*index) ) ) / 2.0 ) ;
        int i = index - 8*j + j*(j-1)/2  + j;
        if (j == i)
            continue;*/
        for(int j = i  + 1; j < 4 ; j++){
            Eigen::VectorXd positive_parameters = parameters;
            positive_parameters(index_map[i]) += h;
            positive_parameters(index_map[j]) += h;
            Eigen::VectorXd negative_parameters = parameters; 
            negative_parameters(index_map[i]) -= h;
            negative_parameters(index_map[j]) -= h;              
            double positive_loglik = calculate_loglikelihood_param(Y_one,  Y_two,  Ones,  lambda,  positive_parameters);
            double negative_loglik = calculate_loglikelihood_param(Y_one,  Y_two,  Ones,  lambda,  negative_parameters);
            hessian(i, j) = hessian(j, i) = (positive_loglik - 2.0*loglik + negative_loglik)/((2.0*h*h) + 0.5*(hessian(i, i) + hessian(j, j)));
        }
       
    }                                                                  
                                                                
} 

static void genetic_correlation_calculate_hessian_and_gradient_multithread(Eigen::MatrixXd & hessian, Eigen::VectorXd & gradient, Eigen::VectorXd Y_one, Eigen::VectorXd Y_two,\
                                                                Eigen::VectorXd Ones,  Eigen::VectorXd lambda, Eigen::VectorXd parameters, double h){
                                                                
double loglik =  calculate_loglikelihood_param(Y_one,  Y_two,  Ones,  lambda,  parameters);

//#pragma omp parallel for
    for(int index = 0; index < 8 ; index++){
        Eigen::VectorXd positive_parameters = parameters;
        positive_parameters(index) += h;
        Eigen::VectorXd negative_parameters = parameters;
        negative_parameters(index) -= h;
        double positive_loglik = calculate_loglikelihood_param(Y_one,  Y_two,  Ones,  lambda,  positive_parameters);
        double negative_loglik = calculate_loglikelihood_param(Y_one,  Y_two,  Ones,  lambda,  negative_parameters);
        gradient(index) =  (positive_loglik - negative_loglik) /(2.0*h);
        hessian(index, index) =  (positive_loglik - 2.0*loglik + negative_loglik) /(h*h);
    }
  
//#pragma omp parallel for
    for(int i = 0 ; i < 8 ; i++){
       /* const long v = 2*8+1;
	    int j = floor( ( v - sqrt( (double)(v*v - 8*index) ) ) / 2.0 ) ;
	    int i = index - 8*j + j*(j-1)/2  + j;
	    if (j == i)
	        continue;*/
        for(int j = i  + 1; j < 8 ; j++){
            Eigen::VectorXd positive_parameters = parameters;
            positive_parameters(i) += h;
            positive_parameters(j) += h;
            Eigen::VectorXd negative_parameters = parameters; 
            negative_parameters(i) -= h;
            negative_parameters(j) -= h;              
            double positive_loglik = calculate_loglikelihood_param(Y_one,  Y_two,  Ones,  lambda,  positive_parameters);
            double negative_loglik = calculate_loglikelihood_param(Y_one,  Y_two,  Ones,  lambda,  negative_parameters);
            hessian(i, j) = hessian(j, i) = (positive_loglik - 2.0*loglik + negative_loglik)/((2.0*h*h) + 0.5*(hessian(i, i) + hessian(j, j)));
        }
       
    }                                                                  
                                                                
}                                                                

static void genetic_correlation_calculate_hessian_and_gradient(Eigen::MatrixXd & hessian, Eigen::VectorXd & gradient, Eigen::VectorXd Y_one, Eigen::VectorXd Y_two,\
                                                                Eigen::VectorXd Ones,  Eigen::VectorXd lambda, double mu_one, double mu_two, double t_one, double t_two,\
                                                                double sd_one, double sd_two, double i_one, double i_two, double h, int constrain_parameter = 0 ){

    double h2r_one = calculate_constraint(t_one);  
    double h2r_two = calculate_constraint(t_two);   
    double rhog = calculate_rho(i_one);  
    double rhoe = calculate_rho(i_two);                                                   
    double mu_one_h = h;
    double mu_two_h = h;
    double t_one_h = h;
    double t_two_h = h;
    double sd_one_h = h;
    double sd_two_h = h;
    double i_one_h = h;
    double i_two_h = h;
  
   double loglik  = calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true);
    gradient(0) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)))/(2.0*mu_one_h) ;
                                                                        
    gradient(1) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two + mu_two_h, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two - mu_two_h, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)))/(2.0*mu_two_h) ;
                                                                 
    gradient(2) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one + t_one_h), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one - t_one_h), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)))/(2.0*t_one_h) ; 
                                                                    
    gradient(3) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two , sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two + t_two_h), calculate_rho(i_one), calculate_rho(i_two)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two - t_two_h), calculate_rho(i_one), calculate_rho(i_two)))/(2.0*t_two_h) ;  
                                                                    
    gradient(4) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two , sd_one + sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one - sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)))/(2.0*sd_one_h) ; 
     
    gradient(5) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two , sd_one, sd_two + sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two - sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two)))/(2.0*sd_two_h) ;    
    
    gradient(6) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two , sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one + i_one_h), calculate_rho(i_two)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one - i_one_h), calculate_rho(i_two)))/(2.0*i_one_h) ;  
 
    gradient(7) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two , sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two + i_two_h)) -calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two - i_two_h)))/(2.0*i_two_h) ;  
   
    
    hessian(0, 0) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two) , true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true))/(mu_one_h*mu_one_h);
    
        hessian(1, 1) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two + h, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two - h, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true))/(mu_two_h*mu_two_h);

    hessian(2, 2) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two, sd_one, sd_two, calculate_constraint(t_one + t_one_h), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one - t_one_h), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true))/(t_one_h*t_one_h);
    
    hessian(3, 3) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two + t_two_h), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two - t_two_h ), calculate_rho(i_one), calculate_rho(i_two), true))/(t_two_h*t_two_h);             

    hessian(4, 4) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two, sd_one + sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one - sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true))/(sd_one_h*sd_one_h);

    hessian(5, 5) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two, sd_one, sd_two + sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two - sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true))/(sd_two_h*sd_two_h);

    hessian(6, 6) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one + i_one_h), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one - i_one_h), calculate_rho(i_two), true))/(i_one_h*i_one_h); 

    hessian(7, 7) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two + i_two_h), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two - i_two_h), true))/(i_two_h*i_two_h);
    
    hessian(0, 1) = hessian(1, 0) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two + mu_two_h, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik+ calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two - mu_two_h, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*mu_one_h*mu_two_h) + 0.5*(hessian(0,0) + hessian(1, 1)));         
        

    hessian(0, 2) = hessian(2, 0) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two , sd_one, sd_two, calculate_constraint(t_one + t_one_h), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik  + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two, calculate_constraint(t_one - t_one_h), calculate_constraint(t_two ), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*t_one_h*mu_one_h) + 0.5*(hessian(0, 0) + hessian(2, 2)));


    hessian(0, 3) = hessian(3, 0) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two , sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two + t_two_h), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two - t_two_h), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*t_two_h*mu_one_h) + 0.5*(hessian(0,0) + hessian(3, 3)));    
    
 

    hessian(0, 4) = hessian(4, 0) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two , sd_one + sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one -sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*sd_one_h*mu_one_h) + 0.5*(hessian(0, 0) + hessian(4, 4)));
    

  
  
      hessian(0, 5) = hessian(5, 0) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two , sd_one, sd_two + sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two -sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*sd_two_h*mu_one_h) + 0.5*(hessian(0, 0) + hessian(5, 5))); 
 

 
     hessian(0, 6) = hessian(6, 0) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two , sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one + i_one_h), calculate_rho(i_two), true) - 2.0*loglik  + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one-  i_one_h), calculate_rho(i_two), true))/((2.0*i_one_h*mu_one_h) + 0.5*(hessian(0, 0) + hessian(6, 6)));  


     hessian(0, 7) = hessian(7, 0) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one + mu_one_h, mu_two , sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two + i_two_h), true) - 2.0*loglik  + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one - mu_one_h, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one ), calculate_rho(i_two - i_two_h), true))/((2.0*i_two_h*mu_one_h) + 0.5*(hessian(0,0) + hessian(7, 7)));     
    


    hessian(1, 2) = hessian(2, 1) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two + mu_two_h, sd_one, sd_two, calculate_constraint(t_one + t_one_h), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik  + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two - mu_two_h, sd_one, sd_two, calculate_constraint(t_one - t_one_h), calculate_constraint(t_two ), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*t_one_h*mu_two_h) + 0.5*(hessian(1, 1) + hessian(2, 2)));       
    


    hessian(1, 3) = hessian(3, 1) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two + mu_two_h, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two + t_two_h), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik  + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two - mu_two_h, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two - t_two_h), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*t_two_h*mu_two_h) + 0.5*(hessian(1, 1) + hessian(3, 3)));   
    

 
    hessian(1, 4) = hessian(4, 1) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two + mu_two_h, sd_one + sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two - mu_two_h, sd_one -sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*sd_one_h*mu_two_h) + 0.5*(hessian(1, 1) + hessian(4, 4))); 
    

    hessian(1, 5) = hessian(5, 1) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two + mu_two_h , sd_one, sd_two + sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik  + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two - mu_two_h, sd_one, sd_two -sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*sd_two_h*mu_two_h) + (hessian(1, 1) + hessian(5, 5)));   
 

 
 
     hessian(1, 6) = hessian(6, 1) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two + mu_two_h, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one + i_one_h), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two - mu_two_h, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one-  i_one_h), calculate_rho(i_two), true))/((2.0*i_one_h*mu_two_h) + 0.5*(hessian(1, 1) + hessian(6, 6)));    


     hessian(1, 7) = hessian(7, 1) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two + mu_two_h , sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two + i_two_h), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two - mu_two_h, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one ), calculate_rho(i_two - i_two_h), true))/((2.0*i_two_h*mu_two_h) + 0.5*(hessian(1, 1) + hessian(7, 7)));   
    
     
    
    hessian(2, 3) = hessian(3, 2) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two, calculate_constraint(t_one + t_one_h), calculate_constraint(t_two + t_two_h), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one - t_one_h), calculate_constraint(t_two - t_two_h), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*t_one_h*t_two_h) + 0.5*(hessian(2, 2) + hessian(3, 3)));
    



    hessian(2, 4) = hessian(4, 2) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one + sd_one_h, sd_two, calculate_constraint(t_one + t_one_h), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik  + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one - sd_one_h, sd_two, calculate_constraint(t_one - t_one_h), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*t_one_h*sd_one_h) + 0.5*(hessian(2, 2) + hessian(4, 4)));
    

  
      hessian(2, 5) = hessian(5, 2) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two + sd_two_h, calculate_constraint(t_one + t_one_h), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two - sd_two_h, calculate_constraint(t_one - t_one_h), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*t_one_h*sd_two_h) + 0.5*(hessian(2, 2) + hessian(5, 5)));    
    

  
      hessian(2, 6) = hessian(6, 2) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two, calculate_constraint(t_one + t_one_h), calculate_constraint(t_two), calculate_rho(i_one +  i_one_h), calculate_rho(i_two), true) - 2.0*loglik  + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one - t_one_h), calculate_constraint(t_two), calculate_rho(i_one - i_one_h), calculate_rho(i_two), true))/((2.0*t_one_h*i_one_h) + 0.5*(hessian(2, 2) + hessian(6, 6)));
    

    hessian(2, 7) = hessian(7, 2) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two, calculate_constraint(t_one + t_one_h), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two + i_two_h), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one - t_one_h), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two - i_two_h), true))/((2.0*t_one_h*i_two_h) + 0.5*(hessian(2, 2) + hessian(7, 7))); 
       
    

    hessian(3, 4) = hessian(4, 3) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one + sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two + t_two_h), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one - sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two - t_two_h), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*t_two_h*sd_one_h) + 0.5*(hessian(3, 3) + hessian(4, 4)));
    


    hessian(3, 5) = hessian(5, 3) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two + sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two + t_two_h), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two - sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two - t_two_h), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*t_two_h*sd_two_h) + 0.5*(hessian(3, 3) + hessian(5, 5)));
    

     hessian(3, 6) = hessian(6, 3) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two + t_two_h), calculate_rho(i_one + i_one_h), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two - t_two_h), calculate_rho(i_one - i_one_h), calculate_rho(i_two), true))/((2.0*t_two_h*i_one_h) + 0.5*(hessian(3, 3) + hessian(6, 6))); 
    


    hessian(3, 7) = hessian(7, 3) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two + t_two_h), calculate_rho(i_one), calculate_rho(i_two + i_two_h), true) - 2.0*loglik  + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two - t_two_h), calculate_rho(i_one), calculate_rho(i_two - i_two_h), true))/((2.0*t_two_h*i_two_h) + 0.5*(hessian(3, 3) + hessian(7, 7)));
    


    hessian(4, 5) = hessian(5, 4) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one + sd_one_h, sd_two + sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one - sd_one_h, sd_two - sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two), true))/((2.0*sd_one_h*sd_two_h) + 0.5*(hessian(4, 4) + hessian(5, 5))); 
    

     hessian(4, 6) = hessian(6, 4) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one + sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two ), calculate_rho(i_one + i_one_h), calculate_rho(i_two), true) - 2.0*loglik  + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one - sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one - i_one_h), calculate_rho(i_two), true))/((2.0*sd_one_h*i_one_h) + 0.5*(hessian(4, 4) + hessian(6, 6))); 
    


    hessian(4, 7) = hessian(7, 4) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one + sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two + i_two_h), true) - 2.0*loglik  + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one - sd_one_h, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two - i_two_h), true))/((2.0*sd_one_h*i_two_h) + 0.5*(hessian(4, 4) + hessian(7,7))); 
    

     hessian(5, 6) = hessian(6, 5) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one, sd_two + sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two ), calculate_rho(i_one + i_one_h), calculate_rho(i_two), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two  - sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one - i_one_h), calculate_rho(i_two), true))/((2.0*sd_two_h*i_one_h) + 0.5*(hessian(5, 5) + hessian(6, 6)));    



     hessian(5, 7) = hessian(7, 5) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one , sd_two + sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two + i_two_h), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two - sd_two_h, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one), calculate_rho(i_two - i_two_h), true))/((2.0*sd_two_h*i_two_h) + 0.5*(hessian(5, 5) + hessian(7, 7)));  


    hessian(6, 7) = hessian(7, 6) = (calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one , mu_two , sd_one , sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one + i_one_h), calculate_rho(i_two + i_two_h), true) - 2.0*loglik + calculate_loglikelihood(Y_one, Y_two, Ones, lambda, mu_one, mu_two, sd_one, sd_two, calculate_constraint(t_one), calculate_constraint(t_two), calculate_rho(i_one - i_one_h), calculate_rho(i_two - i_two_h), true))/((2.0*i_one_h*i_two_h) + 0.5*(hessian(6,6) + hessian(7, 7)));
    

                                                                                                                                           
}


static Eigen::VectorXd calculate_standard_errors(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two,\
                                                                Eigen::MatrixXd covariate_matrix,  Eigen::VectorXd lambda, Eigen::VectorXd parameters, Eigen::VectorXd beta, double h){

Eigen::VectorXd standard_errors;
Eigen::MatrixXd hessian(6 +covariate_matrix.cols()*2, 6 +covariate_matrix.cols()*2);
Eigen::VectorXd all_parameters(6 + covariate_matrix.cols()*2);
//Eigen::VectorXd delta(6 + covariate_matrix.cols()*2);

for(int i = 2; i < 8; i++){
    all_parameters(i - 2) = parameters(i);
    /*
    if(i == 2 || i == 3){
        all_parameters(i - 2) = calculate_constraint(parameters(i));
    }else if (i == 6 || i == 7){
        all_parameters(i - 2) = calculate_rho(parameters(i));
    }else{
        all_parameters(i - 2) = parameters(i);
    }*/
}
for(int i = 0; i < covariate_matrix.cols()*2; i++){
    all_parameters(i + 6) = beta(i);
}
//for(int i = 0; i < all_parameters.rows(); i++){
 //   delta(i) = pow(DP, 0.66667)*max(abs(all_parameters(i)), 1.0);
//}
double loglik =  calculate_loglikelihood_param_two(Y_one,  Y_two,  covariate_matrix,  lambda,  parameters ,beta);
//#pragma omp parallel for
    for(int index = 0; index < (6 + covariate_matrix.cols()*2) ; index++){
      //  Eigen::VectorXd positive_parameters = all_parameters;
      //  positive_parameters(index) += h;
      //  Eigen::VectorXd negative_parameters = all_parameters;
      //  negative_parameters(index) -= h;
      //  double positive_loglik = calculate_loglikelihood_param_all(Y_one,  Y_two,  covariate_matrix,  lambda,  positive_parameters);
      //  double negative_loglik = calculate_loglikelihood_param_all(Y_one,  Y_two,  covariate_matrix,  lambda,  negative_parameters);
       // gradient(index) =  (positive_loglik - negative_loglik) /(2.0*h);
        hessian(index, index) = calculate_hessian(Y_one, Y_two,\
                        covariate_matrix, lambda, \
                        all_parameters,  h, index,\
                         index);// (positive_loglik - 2.0*loglik + negative_loglik) /(h*h);
    }
  
//#pragma omp parallel for
    for(int i = 0 ; i < (6 + covariate_matrix.cols()*2) ; i++){
       /* const long v = 2*8+1;
        int j = floor( ( v - sqrt( (double)(v*v - 8*index) ) ) / 2.0 ) ;
        int i = index - 8*j + j*(j-1)/2  + j;
        if (j == i)
            continue;*/
        for(int j = i  + 1; j < (6 + covariate_matrix.cols()*2) ; j++){
         //   Eigen::VectorXd positive_parameters = all_parameters;
         //   positive_parameters(i) += h;
         //   positive_parameters(j) += h;
          //  Eigen::VectorXd negative_parameters = all_parameters; 
         //   negative_parameters(i) -= h;
         //   negative_parameters(j) -= h;              
         //   double positive_loglik = calculate_loglikelihood_param_all(Y_one,  Y_two,  covariate_matrix,  lambda,  positive_parameters);
         //   double negative_loglik = calculate_loglikelihood_param_all(Y_one,  Y_two,  covariate_matrix,  lambda,  negative_parameters);
            hessian(i, j) = hessian(j, i) = calculate_hessian(Y_one, Y_two,\
                        covariate_matrix, lambda, \
                        all_parameters,  h, i,\
                         j);
            //(positive_loglik - 2.0*loglik + negative_loglik)/((2.0*h*h) + 0.5*(hessian(i, i) + hessian(j, j)));
        }
       
    }
   // cout << "hessian: \n" << hessian << endl;
    standard_errors = (-hessian).inverse().diagonal().cwiseAbs().cwiseSqrt();
    standard_errors(0) =fabs(standard_errors(0)*calculate_dconstraint(all_parameters(0)));
    standard_errors(1) =fabs(standard_errors(1)*calculate_dconstraint(all_parameters(1)));
    standard_errors(4) =fabs(standard_errors(4)*calculate_rho_dconstraint(all_parameters(4)));
    standard_errors(5) =fabs(standard_errors(5)*calculate_rho_dconstraint(all_parameters(5)));    
    return standard_errors;                                                                  
                                                                
} 

/*
static double calculcate_bivariate_model_variate_two(Eigen::VectorXd trait_one, Eigen::VectorXd trait_two, Eigen::VectorXd Ones, Eigen::VectorXd eigenvalues, Eigen::VectorXd & final_parameters, double  & best_h, int & iteration_count, int constrain_parameter = 0, bool debug = false){
    double best_loglikelihood;
    bool first_loglik = false;
    for(double h = 5e-1; h >= 5e-10; h /= 10){
    
        Eigen::VectorXd parameters(7);
        Eigen::VectorXd parameters(8);
        parameters(0) = trait_one.mean();
        parameters(1) = trait_two.mean();
        parameters(2) = 0.35;
        parameters(3) = 0.35;
        parameters(4) = (trait_one - Ones*parameters(0)).norm()/sqrt(Ones.rows());
        parameters(5) = (trait_two - Ones*parameters(1)).norm()/sqrt(Ones.rows());
        parameters(6) = 0.0;
        Eigen::VectorXd delta = Eigen::VectorXd::Zero(7);
       Eigen::MatrixXd hessian(7, 7);
       Eigen::VectorXd gradient(7);        
        int iteration = 0;
       double last_loglik = 0.0;
       double loglik = 0.0;
       double max_delta;              
       do{
        if(debug == true){
            cout << "Differentiation delta: " << h <<  " iteration: " << iteration << endl;
        }
        genetic_correlation_calculate_hessian_and_gradient_version_two( hessian,  gradient,  trait_one,  trait_two,\
                                                                 Ones,  eigenvalues,  parameters(0),  parameters(1),  parameters(2),  parameters(3),\
                                                                 parameters(4),  parameters(5),  parameters(6), h);   
        delta = -hessian.inverse()*gradient;                
        max_delta = -1;
            int param_array[7] = {1, 1, 1, 1, 1, 1, 1};

            double best_loglik = calculate_loglikelihood_two_param(trait_one, trait_two, Ones, eigenvalues, parameters + delta);
            for(int i = 0; i < 7; i++){
                int current_param_array[7] = { 0 };
                current_param_array[i] = 1;
                Eigen::VectorXd current_parameters = parameters;
                current_parameters(i) += delta(i);
               

                
                double current_loglik = calculate_loglikelihood_two_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);
                if(current_loglik > best_loglik){
                    best_loglik = current_loglik;
                             for(int j = 0 ; j < 7 ; j++){
                                param_array[j] = current_param_array[j];
                            }
                }
                
           }
           
           for(int i = 0 ; i < 7; i++){

                for(int param = i+1; param < 7; param++){
                     Eigen::VectorXd current_parameters = parameters;
                 
                    current_parameters(i) += delta(i);
               
              
                    current_parameters(param) += delta(param);
 
                }                                    
                                   
                  //   current_parameters(i) += delta(i);
                   //  current_parameters(param) += delta(param);
                    int current_param_array[7] = { 0 };
                    current_param_array[i] = 1;
                     current_param_array[param] = 1;
                     double current_loglik = calculate_loglikelihood_two_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);
                     
                     if(current_loglik > best_loglik){
                        best_loglik = current_loglik;
                             for(int j = 0 ; j < 7 ; j++){
                                param_array[j] = current_param_array[j];
                            }
                     }
               }
         }
         
         for(int i = 0; i < 7; i++){
            for(int param_i = i+1; param_i < 7 ; param_i++){
                for(int param_j = param_i + 1; param_j < 7 ; param_j++){
                Eigen::VectorXd current_parameters = parameters;   
     
                    current_parameters(i) += delta(i);
              
                    
            
                    current_parameters(param_i) += delta(param_i);
    
 
                    current_parameters(param_j) += delta(param_j);
                       
                    //current_parameters(i) += delta(i);
                   // current_parameters(param_i) += delta(param_i);
                   // current_parameters(param_i) += delta(param_j);     
                    int current_param_array[7] = { 0 };
                    current_param_array[i] = 1;
                    current_param_array[param_i] = 1;  
                    current_param_array[param_j] = 1;
                    double current_loglik = calculate_loglikelihood_two_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);                       
                    if(current_loglik > best_loglik){
                       best_loglik = current_loglik;
                             for(int j = 0 ; j < 7 ; j++){
                                param_array[j] = current_param_array[j];
                            }
                    }
               }
           } 
        }
           
         for(int i = 0; i < 7; i++){
            for(int param_a = i+1; param_a < 7 ; param_a++){
                for(int param_b = param_a + 1; param_b < 7 ; param_b++){
                    for(int param_c = param_b + 1; param_c < 7; param_c++){
                        Eigen::VectorXd current_parameters = parameters;
             
                    current_parameters(i) += delta(i);

  
                    current_parameters(param_a) += delta(param_a);
   
       
                    current_parameters(param_b) += delta(param_b);

     
                    current_parameters(param_c) += delta(param_c);
                                                        
                      //  current_parameters(i) += delta(i);
                     //   current_parameters(param_a) += delta(param_a);
                      //  current_parameters(param_b) += delta(param_b);
                     //   current_parameters(param_c) += delta(param_c);
                        int current_param_array[7] = { 0 };
                        current_param_array[param_a] = 1;
                        current_param_array[param_b] = 1;
                        current_param_array[param_c] = 1;
                        current_param_array[i] = 1;
                        double current_loglik = calculate_loglikelihood_two_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);
                         if(current_loglik > best_loglik){
                             best_loglik = current_loglik;
                             for(int j = 0 ; j < 7 ; j++){
                                param_array[j] = current_param_array[j];
                            }
                         }                        
                        
                    }
               }
           }                                
       }
            
            
       
       
         for(int i = 0; i < 7; i++){
            for(int param_a = i+1; param_a < 7 ; param_a++){
                for(int param_b = param_a + 1; param_b < 7 ; param_b++){
                    for(int param_c = param_b + 1; param_c < 7; param_c++){
                        for(int param_d = param_c + 1; param_d < 7; param_d++){
                        Eigen::VectorXd current_parameters = parameters;
     
                    current_parameters(i) += delta(i);

    
                    current_parameters(param_a) += delta(param_a);
    

                    current_parameters(param_b) += delta(param_b);
 

                    current_parameters(param_c) += delta(param_c);

      
                    current_parameters(param_d) += delta(param_d);
                                          
                       // current_parameters(i) += delta(i);
                      //  current_parameters(param_a) += delta(param_a);
                       // current_parameters(param_b) += delta(param_b);
                      //  current_parameters(param_c) += delta(param_c);
                     //   current_parameters(param_d) += delta(param_d);
                        int current_param_array[7] = { 0 };
                        current_param_array[param_a] = 1;
                        current_param_array[param_b] = 1;
                        current_param_array[param_c] = 1;
                        current_param_array[param_d] = 1;
                        current_param_array[i] = 1;
                        double current_loglik = calculate_loglikelihood_two_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);
                         if(current_loglik > best_loglik){
                             best_loglik = current_loglik;
                             for(int j = 0 ; j < 7 ; j++){
                                param_array[j] = current_param_array[j];
                            }
                         } 
                    }                       
               }
           }                                
         }
        }    
         for(int i = 0; i < 7; i++){
            for(int param_a = i+1; param_a < 7 ; param_a++){
                for(int param_b = param_a + 1; param_b < 7 ; param_b++){
                    for(int param_c = param_b + 1; param_c < 7; param_c++){
                        for(int param_d = param_c + 1; param_d < 7; param_d++){
                            for(int param_e = param_d + 1; param_e < 7; param_e++){
                        Eigen::VectorXd current_parameters = parameters;
                        
        
                    current_parameters(i) += delta(i);

              
                    current_parameters(param_a) += delta(param_a);
   

                    current_parameters(param_b) += delta(param_b);


                    current_parameters(param_c) += delta(param_c);

   
                    current_parameters(param_d) += delta(param_d);
   
  
       
                    current_parameters(param_e) += delta(param_e);

                }                                         
                     //   current_parameters(i) += delta(i);
                    //    current_parameters(param_a) += delta(param_a);
                   //     current_parameters(param_b) += delta(param_b);
                    //    current_parameters(param_c) += delta(param_c);
                   //     current_parameters(param_d) += delta(param_d);
                 //       current_parameters(param_e) += delta(param_e);
                        int current_param_array[7] = { 0 };
                        current_param_array[param_a] = 1;
                        current_param_array[param_b] = 1;
                        current_param_array[param_c] = 1;
                        current_param_array[param_d] = 1;
                        current_param_array[param_e] = 1;
                        current_param_array[i] = 1;
                        double current_loglik = calculate_loglikelihood_two_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);
                         if(current_loglik > best_loglik){
                             best_loglik = current_loglik;
                             for(int j = 0 ; j < 7 ; j++){
                                param_array[j] = current_param_array[j];
                            }
                         } 
                    }                       
               }
           }                                
         }
        }             
       }       
       
         for(int i = 0; i < 7; i++){
            for(int param_a = i+1; param_a < 7 ; param_a++){
                for(int param_b = param_a + 1; param_b < 7 ; param_b++){
                    for(int param_c = param_b + 1; param_c < 7; param_c++){
                        for(int param_d = param_c + 1; param_d < 7; param_d++){
                            for(int param_e = param_d + 1; param_e < 7; param_e++){
                                 for(int param_f = param_e + 1; param_f < 7; param_f++){
                        Eigen::VectorXd current_parameters = parameters;
              
                    current_parameters(i) += delta(i);
            
               
                    current_parameters(param_a) += delta(param_a);
         
    
                    current_parameters(param_b) += delta(param_b);
          
   
                    current_parameters(param_c) += delta(param_c);
              

                    current_parameters(param_d) += delta(param_d);
         
     
                    current_parameters(param_e) += delta(param_e);
  
          
                    current_parameters(param_f) += delta(param_f);
                                        
                        
                     //   current_parameters(i) += delta(i);
                    //    current_parameters(param_a) += delta(param_a);
                    //    current_parameters(param_b) += delta(param_b);
                    //    current_parameters(param_c) += delta(param_c);
                    //    current_parameters(param_d) += delta(param_d);
                    //    current_parameters(param_e) += delta(param_e);
                   //     current_parameters(param_f) += delta(param_f);
                        int current_param_array[7] = { 0 };
                        current_param_array[param_a] = 1;
                        current_param_array[param_b] = 1;
                        current_param_array[param_c] = 1;
                        current_param_array[param_d] = 1;
                        current_param_array[param_e] = 1;
                        current_param_array[param_f] = 1;
                        current_param_array[i] = 1;
                        double current_loglik = calculate_loglikelihood_two_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);
                         if(current_loglik > best_loglik){
                             best_loglik = current_loglik;
                             for(int j = 0 ; j <  ; j++){
                                param_array[j] = current_param_array[j];
                            }
                           
                         } 
                    }                       
               }
           }                                
         }
        }             
       }
      }

              for(int i = 0 ; i < 8; i++){
                if (param_array[i] == 1){
                    if(constrain_parameter == 1 && i == 6)
                        continue;
                    if(constrain_parameter == 2 && i == 7)
                        continue;
                double relative_change = abs(delta(i)/parameters(i)); 
                if(relative_change != relative_change) 
                    continue;            
                parameters(i) += delta(i);
                if ((max_delta < relative_change) || (max_delta == -1 && relative_change != 0))
                       max_delta = relative_change;
            }
                
      }  */ 
static void calculate_mean_and_sd(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two, Eigen::VectorXd Ones, Eigen::VectorXd lambda, Eigen::VectorXd & parameters){
    Eigen::VectorXd omega_one = lambda*calculate_constraint(parameters(2)) + (1.0 - calculate_constraint(parameters(2)))*Eigen::VectorXd::Ones(Y_one.rows());
    omega_one = omega_one.cwiseInverse();
    parameters(0) = Y_one.cwiseProduct(Ones).dot(omega_one)/Ones.cwiseAbs2().dot(omega_one);
    parameters(4) = sqrt((Y_one - Ones*parameters(0)).cwiseAbs2().dot(omega_one)/Y_one.rows());
    Eigen::VectorXd omega_two = lambda*calculate_constraint(parameters(3)) + (1.0 - calculate_constraint(parameters(3)))*Eigen::VectorXd::Ones(Y_one.rows());
    omega_two = omega_two.cwiseInverse();
    parameters(1) = Y_two.cwiseProduct(Ones).dot(omega_two)/Ones.cwiseAbs2().dot(omega_two);
    parameters(5) = sqrt((Y_two - Ones*parameters(1)).cwiseAbs2().dot(omega_two)/Y_two.rows());    

}
static void calculate_mean_and_sd_two(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two, Eigen::MatrixXd covariate_matrix, Eigen::VectorXd lambda, Eigen::VectorXd & parameters, Eigen::VectorXd & beta){
    Eigen::VectorXd omega_one = lambda*calculate_constraint(parameters(2)) + (1.0 - calculate_constraint(parameters(2)))*Eigen::VectorXd::Ones(Y_one.rows());
    Eigen::VectorXd omega_two = lambda*calculate_constraint(parameters(3)) + (1.0 - calculate_constraint(parameters(3)))*Eigen::VectorXd::Ones(Y_one.rows());
    //Eigen::VectorXd omega_one_two = lambda*sqrt(calculate_constraint(parameters(3)) * calculate_constraint(parameters(2)))*calculate_rho(parameters(6)) + \
        Eigen::VectorXd::Ones(lambda.rows())*calculate_rho(parameters(7))*sqrt((1.0-calculate_constraint(parameters(2)))*(1.0-calculate_constraint(parameters(3))));
    omega_one = omega_one.cwiseInverse();
    //Eigen::VectorXd omega = omega_two.cwiseProduct((omega_one.cwiseProduct(omega_two) - omega_one_two.cwiseAbs2()).cwiseInverse());
    Eigen::VectorXd beta_one = (covariate_matrix.transpose()*omega_one.asDiagonal()*covariate_matrix).inverse()*covariate_matrix.transpose()*omega_one.asDiagonal()*Y_one;
    //parameters(0) = Y_one.cwiseProduct(Ones).dot(omega_one)/Ones.cwiseAbs2().dot(omega_one);
    parameters(0) = beta_one(beta_one.rows() - 1);
    parameters(4) = sqrt((Y_one - covariate_matrix*beta_one).cwiseAbs2().dot(omega_one)/Y_one.rows());
    
    omega_two = omega_two.cwiseInverse();
   // omega = omega_one.cwiseProduct((omega_one.cwiseProduct(omega_two) - omega_one_two.cwiseAbs2()).cwiseInverse());
    Eigen::VectorXd beta_two = (covariate_matrix.transpose()*omega_two.asDiagonal()*covariate_matrix).inverse()*covariate_matrix.transpose()*omega_two.asDiagonal()*Y_two;
    parameters(1) = beta_two(beta_two.rows() - 1);
    //parameters(1) = Y_two.cwiseProduct(Ones).dot(omega_two)/Ones.cwiseAbs2().dot(omega_two);
    parameters(5) = sqrt((Y_two - covariate_matrix*beta_two).cwiseAbs2().dot(omega_two)/Y_two.rows());
    for(int i = 0; i < covariate_matrix.cols(); i++){
        beta(i) = beta_one(i);
        beta(i + covariate_matrix.cols()) = beta_two(i);
    }    

}
static void calculate_beta(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two, Eigen::MatrixXd covariate_matrix, Eigen::VectorXd lambda, Eigen::VectorXd & parameters, Eigen::VectorXd & beta){
    const double h2_one = calculate_constraint(parameters(2));
    const double h2_two = calculate_constraint(parameters(3));
    const double sd_one = abs(parameters(4));
    const double sd_two = abs(parameters(5));
    const double rhog = calculate_rho(parameters(6));
    const double rhoe = calculate_rho(parameters(7));
    Eigen::VectorXd omega_one_two = (lambda*sqrt(h2_one*h2_two)*rhog + Eigen::VectorXd::Ones(lambda.rows())*sqrt((1.0-h2_two)*(1.0-h2_one))*rhoe)*sd_one*sd_two;         
    Eigen::VectorXd omega_one = (lambda*h2_one + (1.0 - h2_one)*Eigen::VectorXd::Ones(Y_one.rows()))*pow(sd_one, 2);
    Eigen::VectorXd omega_two = (lambda*h2_two + (1.0 - h2_two)*Eigen::VectorXd::Ones(Y_one.rows()))*pow(sd_two, 2);
    compute_blockwise_diagonal_inversions(omega_one, omega_one_two, omega_two);
    Eigen::MatrixXd m = covariate_matrix.transpose()*omega_one_two.asDiagonal()*covariate_matrix;
    Eigen::MatrixXd v = covariate_matrix.transpose()*omega_two.asDiagonal()*covariate_matrix;
    Eigen::MatrixXd n = covariate_matrix.transpose()*omega_one.asDiagonal()*covariate_matrix;
    Eigen::VectorXd w = covariate_matrix.transpose()*omega_one.asDiagonal()*Y_one;
    Eigen::VectorXd o = covariate_matrix.transpose()*omega_one_two.asDiagonal()*Y_one;
    Eigen::VectorXd h = covariate_matrix.transpose()*omega_one_two.asDiagonal()*Y_two;
    Eigen::VectorXd q = covariate_matrix.transpose()*omega_two.asDiagonal()*Y_two;
    Eigen::VectorXd beta_one = (-m*m + n*v).inverse()*(v*w +  v*h + m*(-q-o));
    Eigen::VectorXd beta_two = -(-m*m + n*v).inverse()*(m*w + n*(-q-o) + m*h);
    parameters(0) = beta_one(beta_one.rows() - 1);
    
    
   
    parameters(1) = beta_two(beta_two.rows() - 1);

    for(int i = 0; i < covariate_matrix.cols(); i++){
        beta(i) = beta_one(i);
        beta(i + covariate_matrix.cols()) = beta_two(i);

    }                  
}
static void calculate_mean_and_sd_three(Eigen::VectorXd Y_one, Eigen::VectorXd Y_two, Eigen::MatrixXd covariate_matrix, Eigen::VectorXd lambda, Eigen::VectorXd & parameters, Eigen::VectorXd & beta){
    double rhop =parameters(6)*sqrt(parameters(2))*\
             sqrt(parameters(3)) +parameters(7)*\
             sqrt(1.0-parameters(2))*sqrt(1.0-parameters(3));

    Eigen::VectorXd omega_one_two = (lambda*sqrt(parameters(2))*sqrt(parameters(3))*parameters(6) + Eigen::VectorXd::Ones(lambda.rows())*sqrt(1.0-parameters(2))*sqrt(1.0-parameters(3))*parameters(7));//parameters(4)*parameters(5);         
    Eigen::VectorXd omega_one = (lambda*parameters(2) + (1.0 - parameters(2))*Eigen::VectorXd::Ones(Y_one.rows()));//pow(parameters(4), 2);
    Eigen::VectorXd omega_two = (lambda*parameters(3) + (1.0 - parameters(3))*Eigen::VectorXd::Ones(Y_one.rows()));//pow(parameters(5), 2);
    //compute_blockwise_diagonal_inversions(omega_one, omega_one_two, omega_two);
    omega_one = omega_one.cwiseInverse();
    omega_two = omega_two.cwiseInverse();
   //omega_one = omega_one.cwiseInverse()*pow(1-rhop*rhop, -1);
    // omega_two = omega_two.cwiseInverse()*pow(1-rhop*rhop, -1);
    //omega_one_two = omega_one_two.cwiseInverse()*pow(1-rhop*rhop, -1);
   //  -(x*s^-1*x)^-1*(((x*k^-1*p) - (x*k^-1*x*((x*r^-1*x)^-1)*(-((x*k^-1*y) - (x*k^-1*x*b)) + (x*r^-1*p)))) + (x*s^-1*y)))= (x*s^-1*x*b)) 
   //  Eigen::VectorXd beta_one=  ((covariate_matrix.transpose()*omega_one.asDiagonal()*covariate_matrix).inverse()* \
        (covariate_matrix.transpose()*omega_one.asDiagonal()*Y_one -((covariate_matrix.transpose()*omega_one_two.asDiagonal()*Y_two) - (covariate_matrix.transpose()*omega_one_two.asDiaognal()*covariate_matrix*((covariate_matrix.transpose()*omega_two.asDiagonal()*covariate_matrix).inverse()*\
            (-(covariate_matrix.transpose()*omega_one_two*Y_one) - (covariate_matrix.transpose()*omege_one_two.asDiagonal()*)
   /* Eigen::MatrixXd m = covariate_matrix.transpose()*omega_one_two.asDiagonal()*covariate_matrix;
    Eigen::MatrixXd v = covariate_matrix.transpose()*omega_two.asDiagonal()*covariate_matrix;
    Eigen::MatrixXd n = covariate_matrix.transpose()*omega_one.asDiagonal()*covariate_matrix;
    Eigen::VectorXd w = covariate_matrix.transpose()*omega_one.asDiagonal()*Y_one;
    Eigen::VectorXd o = covariate_matrix.transpose()*omega_one_two.asDiagonal()*Y_one;
    Eigen::VectorXd h = covariate_matrix.transpose()*omega_one_two.asDiagonal()*Y_two;
    Eigen::VectorXd q = covariate_matrix.transpose()*omega_two.asDiagonal()*Y_two;*/
   // Eigen::VectorXd beta_one = (-m*m + n*v).inverse()*(v*w +  v*h + m*(-q-o));
  //  Eigen::VectorXd beta_one = (n-m*m*v.inverse()).inverse()*(m*v.inverse()*(q-o) - h + w);
   // Eigen::VectorXd beta_two = -(-m*m + n*v).inverse()*(m*w + n*(-q-o) + m*h);
    //Eigen::VectorXd beta_one = (n*v-m*m).inverse()*(v*w - v*h + m*(q-o));
   // Eigen::VectorXd beta_two = (n*v - m*m).inverse()*(m*w + n*(q-o) - m*h); 
   // Eigen::VectorXd beta_one = (covariate_matrix.transpose()*omega_one.asDiagonal()*covariate_matrix).inverse()*covariate_matrix.transpose()*omega_one.asDiagonal()*Y_one;
    //parameters(0) = Y_one.cwiseProduct(Ones).dot(omega_one)/Ones.cwiseAbs2().dot(omega_one);
    Eigen::VectorXd beta_one = (covariate_matrix.transpose()*omega_one.asDiagonal()*covariate_matrix).inverse()*covariate_matrix.transpose()*omega_one.asDiagonal()*Y_one;
    Eigen::VectorXd beta_two = (covariate_matrix.transpose()*omega_two.asDiagonal()*covariate_matrix).inverse()*covariate_matrix.transpose()*omega_two.asDiagonal()*Y_two;
    parameters(0) = beta_one(beta_one.rows() - 1);
    parameters(4) = sqrt((Y_one - covariate_matrix*beta_one).cwiseAbs2().dot(omega_one)/Y_one.rows());
    
   
    parameters(1) = beta_two(beta_two.rows() - 1);
    //parameters(1) = Y_two.cwiseProduct(Ones).dot(omega_two)/Ones.cwiseAbs2().dot(omega_two);
    parameters(5) = sqrt((Y_two - covariate_matrix*beta_two).cwiseAbs2().dot(omega_two)/Y_two.rows());
    for(int i = 0; i < covariate_matrix.cols(); i++){
        beta(i) = beta_one(i);
        beta(i + covariate_matrix.cols()) = beta_two(i);

    }    

}

static double calculate_bivariate_model_three(Eigen::VectorXd trait_one, Eigen::VectorXd trait_two, Eigen::MatrixXd covariate_matrix, \
    Eigen::VectorXd eigenvalues, Eigen::VectorXd & final_parameters, Eigen::VectorXd & final_beta, Eigen::VectorXd & standard_errors, \
    double  & best_h, int & iteration_count, int constrain_parameter = 0, bool debug = false){
/*
    Eigen::VectorXd fixed_parameters(8);
    
    fixed_parameters(0) = 1895.377850175171;
    fixed_parameters(1) = 1807.993098790394;
    fixed_parameters(2) = 0.09675017378418029;
    fixed_parameters(3) = 0.2603526383507959;
    fixed_parameters(4) = 483.5970535393554;
    fixed_parameters(5) = 514.1873061645562;
    fixed_parameters(6) = 0.9341313728734818;
    fixed_parameters(7) = 0.1677193557134624;
    fixed_parameters(0) = 1895.757526045272;
    fixed_parameters(1) = 1807.993098790394;
    fixed_parameters(2) = 0.1028159340403367;
    fixed_parameters(3) = 0.261393096986652;
    fixed_parameters(4) = 484.5567230587188;
    fixed_parameters(5) = 517.2392210735813;
    fixed_parameters(6) = 0.9371208575789899;
    fixed_parameters(7) = 0.1674260093316723;   
    Eigen::VectorXd fixed_beta(covariate_matrix.cols()*2);
    fixed_beta(0) = -8.386926261594702;
    fixed_beta(2) = -15.59735530025079;
    fixed_beta(1) = 1895.757526045272;
    fixed_beta(3) = 1807.993098790394;
    double first_fixed_loglik = calculate_loglikelihood_param_two(trait_one, trait_two, covariate_matrix, eigenvalues, fixed_parameters, fixed_beta, false);
    
    calculate_mean_and_sd_three(trait_one, trait_two, covariate_matrix, eigenvalues, fixed_parameters , fixed_beta);
    double second_fixed_loglik = calculate_loglikelihood_param_two(trait_one, trait_two, covariate_matrix, eigenvalues, fixed_parameters, fixed_beta, false);
    cout << first_fixed_loglik << endl;
    cout << second_fixed_loglik << endl;
    cout << "hello\n";
    cout << "beta: \n";
    cout << fixed_beta << endl;
    cout << "sd: \n";
    cout << fixed_parameters(4) << " " << fixed_parameters(5) << endl;
    return 0.0;*/
    double best_loglikelihood;
    bool first_loglik = false;
    double h;
  //  omp_set_num_threads(5);
    const int index_map[4] = {2, 3, 6, 7};
//#pragma omp parallel for private(h) shared(first_loglik)
    for(int i = 2; i <= 2 ; i++){
        if (i == 1){ 
         cout << "Using " << omp_get_num_threads() << " threads for computation\n"; 
        }
        h = best_h;//pow(10, -i);
        
  //  for(double h = 5e-1; h >= 5e-10; h /= 10){
       Eigen::VectorXd beta(covariate_matrix.cols()*2);
       Eigen::VectorXd parameters(8);
       parameters(2) = 0.8;
       parameters(3) = 0.8;
       parameters(6) = 0.0;
       parameters(7) = 0.0;
       calculate_mean_and_sd_two(trait_one, trait_two, covariate_matrix, eigenvalues, parameters ,beta);
       double loglik = calculate_loglikelihood_param_two(trait_one, trait_two, covariate_matrix, eigenvalues, parameters, beta);
    /*   parameters(0) = trait_one.mean();
       parameters(1) = trait_two.mean();
       parameters(2) = 0.35;
       parameters(3) = 0.35;
       parameters(4) = (trait_one - Ones*parameters(0)).norm()/sqrt(Ones.rows());
       parameters(5) = (trait_two - Ones*parameters(1)).norm()/sqrt(Ones.rows());
       parameters(6) = 0.0;
       parameters(7) = 0.0;*/
       Eigen::VectorXd delta = Eigen::VectorXd::Zero(8);
       Eigen::MatrixXd hessian(4, 4);
       Eigen::VectorXd gradient(4);
       int iteration = 0;
       double last_loglik = 0.0;
       double loglik_error;
       //double loglik = 0.0;
       double max_delta;
       genetic_correlation_calculate_hessian_and_gradient_multithread_three(hessian,  gradient, trait_one, trait_two,\
                                                               covariate_matrix,  eigenvalues,  parameters, beta, h);
       delta = -hessian.inverse()*gradient;
       double error;
       int converge_count = 0;                                                               
       do{
        if(debug == true){
            cout << "Differentiation delta: " << h <<  " iteration: " << iteration << endl;
        }
       // genetic_correlation_calculate_hessian_and_gradient( hessian,  gradient,  trait_one,  trait_two,\
                                                                 Ones,  eigenvalues,  parameters(0),  parameters(1),  parameters(2),  parameters(3),\
                                                                 parameters(4),  parameters(5),  parameters(6),  parameters(7),  h, constrain_parameter);  

                                                                                                             
       // delta = -hessian.inverse()*gradient;                
        //max_delta = -1;
            int param_array[4] = {1, 1, 1, 1};
            if(constrain_parameter == 1){
                delta(2) = 0.0;
            }
            if(constrain_parameter == 2){
                delta(3) = 0.0;
            }
            Eigen::VectorXd test_parameters = parameters;
            test_parameters(2) += delta(0);
            test_parameters(3) += delta(1);
            test_parameters(6) += delta(2);
            test_parameters(7) += delta(3);
            Eigen::VectorXd test_beta = beta;
            calculate_mean_and_sd_two(trait_one, trait_two, covariate_matrix, eigenvalues, test_parameters, test_beta);
            double best_loglik = calculate_loglikelihood_param_two(trait_one, trait_two, covariate_matrix, eigenvalues, test_parameters, test_beta);
            Eigen::VectorXd best_parameters = test_parameters;
            Eigen::VectorXd best_beta = test_beta;
            for(int step = 1; step <= 19 ; step++){
               Eigen::VectorXd current_parameters = parameters;
               Eigen::VectorXd current_beta = beta;
               current_parameters(2) += step*0.05*delta(0);
               current_parameters(3) += step*0.05*delta(1);
               current_parameters(6) += step*0.05*delta(2);
               current_parameters(7) += step*0.05*delta(3);
               calculate_mean_and_sd_two(trait_one, trait_two, covariate_matrix, eigenvalues, current_parameters, current_beta);
               double current_loglik = calculate_loglikelihood_param_two(trait_one, trait_two, covariate_matrix, eigenvalues, current_parameters, current_beta);
               if(current_loglik > best_loglik || best_loglik == 0.0){
                   best_loglik = current_loglik;
                   best_beta = current_beta;
                   best_parameters = current_parameters;
                   //for(int j = 0; j < 4; j++){
                    //   param_array[j] = current_param_array[j];
                   //}
               }               

            }
            /*for(unsigned param_index = 1; param_index < 15 ; param_index++){
               int current_param_array[4] = {0};
               Eigen::VectorXd current_parameters = parameters;
               Eigen::VectorXd current_beta = beta;
               if((param_index & 1) == 0){
                   current_parameters(2) += delta(0);
                   current_param_array[0] = 1;
               }
               if((param_index & 2) == 0){
                   current_parameters(3) += delta(1);
                   current_param_array[1] = 1;
               }
               if((param_index & 4) == 0){
                   current_parameters(6) += delta(2);
                   current_param_array[2] = 1;
               }
               if((param_index & 8) == 0){
                   current_parameters(7) += delta(3);
                   current_param_array[3] = 1;
               }
               calculate_mean_and_sd_two(trait_one, trait_two, covariate_matrix, eigenvalues, current_parameters, current_beta);
               double current_loglik = calculate_loglikelihood_param_two(trait_one, trait_two, covariate_matrix, eigenvalues, current_parameters, current_beta);
               if(current_loglik > best_loglik || best_loglik == 0.0){
                   best_loglik = current_loglik;
                   best_beta = current_beta;
                   best_parameters = current_parameters;
                   for(int j = 0; j < 4; j++){
                       param_array[j] = current_param_array[j];
                   }
               }

               
           }*/
           beta = best_beta;
           last_loglik = loglik;
           loglik = best_loglik;
           loglik_error = fabs((last_loglik - loglik)/last_loglik);
           parameters = best_parameters;
           /*
              for(int i = 0 ; i < 4; i++){
                if (param_array[i] == 1){
                    if(constrain_parameter == 1 && i == 2)
                        continue;
                    if(constrain_parameter == 2 && i == 3)
                        continue;
                double relative_change = 0;// = abs(delta(i));
                if (delta(i) != delta(i))
                    continue; 
               // if(relative_change != relative_change) 
               //     continue;
                if(i == 0 || i == 1){
                    double old_parameter = parameters(index_map[i]);
                    double old_h2r = calculate_constraint(old_parameter);
                    parameters(index_map[i]) += delta(i);
                    double new_h2r = calculate_constraint(parameters(index_map[i]));
                    //relative_change = abs(new_h2r - old_h2r);

                }else{
                    double old_parameter = parameters(index_map[i]);
                    double old_rho = calculate_rho(old_parameter);
                    parameters(index_map[i]) += delta(i);
                    double new_rho = calculate_rho(parameters(index_map[i]));
                    //relative_change = abs(new_rho - old_rho);                    
                }                
               // parameters(index_map[i]) += delta(i);
                if ((max_delta < relative_change) || (max_delta == -1 && relative_change != 0))
                       max_delta = relative_change;
            }
                
      }  
      
     calculate_mean_and_sd_two(trait_one, trait_two, covariate_matrix, eigenvalues, parameters, beta);

      
         last_loglik = loglik; 
         loglik =  calculate_loglikelihood_param_two(trait_one, trait_two, covariate_matrix, eigenvalues, parameters, beta);                                                            

     */                                                                                                                                                                                                                                                                                                   

        genetic_correlation_calculate_hessian_and_gradient_multithread_three(hessian,  gradient, trait_one, trait_two,\
                                                               covariate_matrix,  eigenvalues,  parameters, beta, h);
        delta = -hessian.inverse()*gradient;

        error = fabs(delta.dot(gradient));
      if (debug == true){
        double h2r_one = calculate_constraint(parameters(2)); 
        double h2r_two = calculate_constraint(parameters(3)); 
        double rhog = calculate_rho(parameters(6)); 
        double rhoe = calculate_rho(parameters(7));        
        cout << "iteration: " << iteration << endl;
        cout << "mean one: " << parameters(0) << " mean two: " << parameters(1) << endl;
        cout << "sd one: " << parameters(4) << " sd two: " << parameters(5) << endl;
        cout << "t one: " << parameters(2) << " h2r one: " << h2r_one << endl;
        cout << "t two: " << parameters(3) << " h2r two: " << h2r_two << endl;
        cout << "i one: " << parameters(6) << " rhog: " << rhog << endl;
        cout << "i two: " << parameters(7) << " rhoe: " << rhoe << endl;
        cout << "loglik: " << loglik << endl;

     }
        if (loglik_error <= MAX_LOGLIK_ERROR){
            converge_count++;
        }else{
            converge_count = 0;
        } 
        /*
        if(converge_count == 3 && h != 1E-06){
            converge_count = 0;
            h /= 10;
        }  */     
        iteration++; 
    }while(loglik == loglik && delta == delta && iteration < MAX_ITERATIONS  && (converge_count  < 3));
    #pragma omp critical 
    {
    if(loglik == loglik && delta == delta && (loglik > best_loglikelihood || first_loglik == false)){

        best_h = h;
        best_loglikelihood = loglik;
        final_parameters = parameters;
        final_beta = beta;
        first_loglik = true;
        iteration_count = iteration;
        standard_errors = calculate_standard_errors(trait_one, trait_two,\
                                        covariate_matrix,  eigenvalues, final_parameters, final_beta,  best_h);
     }
    }
  }
  if(first_loglik == false){
    standard_errors = Eigen::VectorXd::Zero(6+covariate_matrix.cols()*2);
    return nan("");
  }

  
return best_loglikelihood;

  
} 
static double calculate_bivariate_model_two(Eigen::VectorXd trait_one, Eigen::VectorXd trait_two, Eigen::VectorXd Ones, Eigen::VectorXd eigenvalues, Eigen::VectorXd & final_parameters, double  & best_h, int & iteration_count, int constrain_parameter = 0, bool debug = false){

    double best_loglikelihood;
    bool first_loglik = false;
    double h;
    omp_set_num_threads(10);
    const int index_map[4] = {2, 3, 6, 7};
#pragma omp parallel for private(h) shared(first_loglik)
    for(int i = 2; i <= 11 ; i++){
        if (i == 1){ 
         cout << "Using " << omp_get_num_threads() << " threads for computation\n"; 
        }
        h = pow(10, -i);
        
  //  for(double h = 5e-1; h >= 5e-10; h /= 10){
      
       Eigen::VectorXd parameters(8);
       parameters(2) = 0.35;
       parameters(3) = 0.35;
       parameters(6) = 0.0;
       parameters(7) = 0.0;
       calculate_mean_and_sd(trait_one, trait_two, Ones, eigenvalues, parameters);
    /*   parameters(0) = trait_one.mean();
       parameters(1) = trait_two.mean();
       parameters(2) = 0.35;
       parameters(3) = 0.35;
       parameters(4) = (trait_one - Ones*parameters(0)).norm()/sqrt(Ones.rows());
       parameters(5) = (trait_two - Ones*parameters(1)).norm()/sqrt(Ones.rows());
       parameters(6) = 0.0;
       parameters(7) = 0.0;*/
       Eigen::VectorXd delta = Eigen::VectorXd::Zero(8);
       Eigen::MatrixXd hessian(4, 4);
       Eigen::VectorXd gradient(4);
       int iteration = 0;
       double last_loglik = 0.0;
       double loglik = 0.0;
       double max_delta;
       do{
        if(debug == true){
            cout << "Differentiation delta: " << h <<  " iteration: " << iteration << endl;
        }
       // genetic_correlation_calculate_hessian_and_gradient( hessian,  gradient,  trait_one,  trait_two,\
                                                                 Ones,  eigenvalues,  parameters(0),  parameters(1),  parameters(2),  parameters(3),\
                                                                 parameters(4),  parameters(5),  parameters(6),  parameters(7),  h, constrain_parameter);  
        genetic_correlation_calculate_hessian_and_gradient_multithread_two(hessian,  gradient, trait_one, trait_two,\
                                                               Ones,  eigenvalues,  parameters,  h);
                                                                                                             
        delta = -hessian.inverse()*gradient;                
        max_delta = -1;
            int param_array[4] = {1, 1, 1, 1};
            if(constrain_parameter == 1){
                delta(2) = 0.0;
            }
            if(constrain_parameter == 2){
                delta(3) = 0.0;
            }
            Eigen::VectorXd test_parameters = parameters;
            test_parameters(2) += delta(0);
            test_parameters(3) += delta(1);
            test_parameters(6) += delta(2);
            test_parameters(7) += delta(3);
            calculate_mean_and_sd(trait_one, trait_two, Ones, eigenvalues, test_parameters);
            double best_loglik = calculate_loglikelihood_param(trait_one, trait_two, Ones, eigenvalues, test_parameters);

            for(unsigned param_index = 1; param_index < 15 ; param_index++){
               int current_param_array[4] = {0};
               Eigen::VectorXd current_parameters = parameters;
               if((param_index & 1) == 0){
                   current_parameters(2) += delta(0);
                   current_param_array[0] = 1;
               }
               if((param_index & 2) == 0){
                   current_parameters(3) += delta(1);
                   current_param_array[1] = 1;
               }
               if((param_index & 4) == 0){
                   current_parameters(6) += delta(2);
                   current_param_array[2] = 1;
               }
               if((param_index & 8) == 0){
                   current_parameters(7) += delta(3);
                   current_param_array[3] = 1;
               }
               calculate_mean_and_sd(trait_one, trait_two, Ones, eigenvalues, current_parameters);
               double current_loglik = calculate_loglikelihood_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);
               if(current_loglik > best_loglik || best_loglik == 0.0){
                   best_loglik = current_loglik;
                   for(int j = 0; j < 4; j++){
                       param_array[j] = current_param_array[j];
                   }
               }

               
           }


              for(int i = 0 ; i < 4; i++){
                if (param_array[i] == 1){
                    if(constrain_parameter == 1 && i == 2)
                        continue;
                    if(constrain_parameter == 2 && i == 3)
                        continue;
                double relative_change = 0;// = abs(delta(i));
                if (delta(i) != delta(i))
                    continue; 
               // if(relative_change != relative_change) 
               //     continue;
                if(i == 0 || i == 1){
                    double old_parameter = parameters(index_map[i]);
                    double old_h2r = calculate_constraint(old_parameter);
                    parameters(index_map[i]) += delta(i);
                    double new_h2r = calculate_constraint(parameters(index_map[i]));
                    relative_change = abs(new_h2r - old_h2r);

                }else{
                    double old_parameter = parameters(index_map[i]);
                    double old_rho = calculate_rho(old_parameter);
                    parameters(index_map[i]) += delta(i);
                    double new_rho = calculate_rho(parameters(index_map[i]));
                    relative_change = abs(new_rho - old_rho);                    
                }                
               // parameters(index_map[i]) += delta(i);
                if ((max_delta < relative_change) || (max_delta == -1 && relative_change != 0))
                       max_delta = relative_change;
            }
                
      }  
      
     calculate_mean_and_sd(trait_one, trait_two, Ones, eigenvalues, parameters);
     double h2r_one = calculate_constraint(parameters(2)); 
     double h2r_two = calculate_constraint(parameters(3)); 
     double rhog = calculate_rho(parameters(6)); 
     double rhoe = calculate_rho(parameters(7));
      
         last_loglik = loglik; 
         loglik =  calculate_loglikelihood_param(trait_one, trait_two, Ones, eigenvalues, parameters);                                                            

                                                                                                                                                                                                                                                                                                        
      if (debug == true){

        cout << "mean one: " << parameters(0) << " mean two: " << parameters(1) << endl;
        cout << "sd one: " << parameters(4) << " sd two: " << parameters(5) << endl;
        cout << "t one: " << parameters(2) << " h2r one: " << h2r_one << endl;
        cout << "t two: " << parameters(3) << " h2r two: " << h2r_two << endl;
        cout << "i one: " << parameters(6) << " rhog: " << rhog << endl;
        cout << "i two: " << parameters(7) << " rhoe: " << rhoe << endl;
        cout << "loglik: " << loglik << endl;

     } 

        iteration++; 
    }while(loglik == loglik && delta == delta && iteration < MAX_ITERATIONS  && (max_delta >= MAX_DELTA_ERROR));
    #pragma omp critical 
    {
    if(loglik == loglik && delta == delta && (loglik > best_loglikelihood || first_loglik == false)){
        best_h = h;
        best_loglikelihood = loglik;
        final_parameters = parameters;
        first_loglik = true;
        iteration_count = iteration;
     }
    }
  }
  if(first_loglik == false){
    return nan("");
  }
  return best_loglikelihood;

  
}           
static double calculate_bivariate_model(Eigen::VectorXd trait_one, Eigen::VectorXd trait_two, Eigen::VectorXd Ones, Eigen::VectorXd eigenvalues, Eigen::VectorXd & final_parameters, double  & best_h, int & iteration_count, int constrain_parameter = 0, bool debug = false){

    double best_loglikelihood;
    bool first_loglik = false;
    double h;
    omp_set_num_threads(10);
   
#pragma omp parallel for private(h) shared(first_loglik)
    for(int i = 1; i <= 10 ; i++){
        if (i == 1){ 
         cout << "Using " << omp_get_num_threads() << " threads for computation\n"; 
        }
        h = pow(10, -i/2.0);
        
  //  for(double h = 5e-1; h >= 5e-10; h /= 10){
      
       Eigen::VectorXd parameters(8);
       parameters(0) = trait_one.mean();
       parameters(1) = trait_two.mean();
       parameters(2) = 0.35;
       parameters(3) = 0.35;
       parameters(4) = (trait_one - Ones*parameters(0)).norm()/sqrt(Ones.rows());
       parameters(5) = (trait_two - Ones*parameters(1)).norm()/sqrt(Ones.rows());
       parameters(6) = 0.0;
       parameters(7) = 0.0;
       Eigen::VectorXd delta = Eigen::VectorXd::Zero(8);
       Eigen::MatrixXd hessian(8, 8);
       Eigen::VectorXd gradient(8);
       int iteration = 0;
       double last_loglik = 0.0;
       double loglik = 0.0;
       double max_delta;
       do{
        if(debug == true){
            cout << "Differentiation delta: " << h <<  " iteration: " << iteration << endl;
        }
       // genetic_correlation_calculate_hessian_and_gradient( hessian,  gradient,  trait_one,  trait_two,\
                                                                 Ones,  eigenvalues,  parameters(0),  parameters(1),  parameters(2),  parameters(3),\
                                                                 parameters(4),  parameters(5),  parameters(6),  parameters(7),  h, constrain_parameter);  
        genetic_correlation_calculate_hessian_and_gradient_multithread(hessian,  gradient, trait_one, trait_two,\
                                                               Ones,  eigenvalues,  parameters,  h);                                                         
        delta = -hessian.inverse()*gradient;                
        max_delta = -1;
            int param_array[8] = {1, 1, 1, 1, 1, 1, 1, 1};
            if(constrain_parameter == 1){
                delta(6) = 0.0;
            }
            if(constrain_parameter == 2){
                delta(7) = 0.0;
            }
            double best_loglik = calculate_loglikelihood_param(trait_one, trait_two, Ones, eigenvalues, parameters + delta);
            for(unsigned param_index = 0; param_index < 256 ; param_index++){
               int current_param_array[8] = {0};
               Eigen::VectorXd current_parameters = parameters;
               if((param_index & 1) == 0){
                   current_parameters(0) += delta(0);
                   current_param_array[0] = 1;
               }
               if((param_index & 2) == 0){
                   current_parameters(1) += delta(1);
                   current_param_array[1] = 1;
               }
               if((param_index & 4) == 0){
                   current_parameters(2) += delta(2);
                   current_param_array[2] = 1;
               }
               if((param_index & 8) == 0){
                   current_parameters(3) += delta(3);
                   current_param_array[3] = 1;
               }
               if((param_index & 16) == 0){
                   current_parameters(4) += delta(4);
                   current_param_array[4] = 1;
               }
               if((param_index & 32) == 0){
                   current_parameters(5) += delta(5);
                   current_param_array[5] = 1;
               }
               if((param_index & 64) == 0 ){
                   current_parameters(6) += delta(6);
                   current_param_array[6] = 1;
               }
               if((param_index & 128) == 0 ){
                   current_parameters(7) += delta(7);
                   current_param_array[7] = 1;
               }
               double current_loglik = calculate_loglikelihood_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);
               if(current_loglik > best_loglik || best_loglik == 0.0){
                   best_loglik = current_loglik;
                   for(int j = 0; j < 8; j++){
                       param_array[j] = current_param_array[j];
                   }
               }

               
           }
/*
            for(int i = 0; i < 8; i++){
                int current_param_array[8] = { 0 };
                current_param_array[i] = 1;
                Eigen::VectorXd current_parameters = parameters;
                if (constrain_parameter == 0){
                    current_parameters(i) += delta(i);
                }else if (constrain_parameter == 1 && i != 6){
                   current_parameters(i) += delta(i);
                }else if (constrain_parameter == 2 && i != 7){ 
                   current_parameters(i) += delta(i); 
                }               

                
                double current_loglik = calculate_loglikelihood_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);
                if(current_loglik > best_loglik){
                    best_loglik = current_loglik;
                             for(int j = 0 ; j < 8 ; j++){
                                param_array[j] = current_param_array[j];
                            }
                }
                
           }
           
           for(int i = 0 ; i < 8; i++){

                for(int param = i+1; param < 8; param++){
                     Eigen::VectorXd current_parameters = parameters;
                 if (constrain_parameter == 0){
                    current_parameters(i) += delta(i);
                }else if (constrain_parameter == 1 && i != 6){
                   current_parameters(i) += delta(i);
                }else if (constrain_parameter == 2 && i != 7){ 
                   current_parameters(i) += delta(i); 
                }
                if (constrain_parameter == 0){
                    current_parameters(param) += delta(param);
                }else if (constrain_parameter == 1 && param != 6){
                   current_parameters(param) += delta(param);
                }else if (constrain_parameter == 2 && param != 7){ 
                   current_parameters(param) += delta(param); 
                }                                    
                                   
                  //   current_parameters(i) += delta(i);
                   //  current_parameters(param) += delta(param);
                    int current_param_array[8] = { 0 };
                    current_param_array[i] = 1;
                     current_param_array[param] = 1;
                     double current_loglik = calculate_loglikelihood_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);
                     
                     if(current_loglik > best_loglik){
                        best_loglik = current_loglik;
                             for(int j = 0 ; j < 8 ; j++){
                                param_array[j] = current_param_array[j];
                            }
                     }
               }
         }
         
         for(int i = 0; i < 8; i++){
            for(int param_i = i+1; param_i < 8 ; param_i++){
                for(int param_j = param_i + 1; param_j < 8 ; param_j++){
                Eigen::VectorXd current_parameters = parameters;   
                 if (constrain_parameter == 0){
                    current_parameters(i) += delta(i);
                }else if (constrain_parameter == 1 && i != 6){
                   current_parameters(i) += delta(i);
                }else if (constrain_parameter == 2 && i != 7){ 
                   current_parameters(i) += delta(i); 
                }                
                    
                 if (constrain_parameter == 0){
                    current_parameters(param_i) += delta(param_i);
                }else if (constrain_parameter == 1 && param_i != 6){
                   current_parameters(param_i) += delta(param_i);
                }else if (constrain_parameter == 2 && param_i != 7){ 
                   current_parameters(param_i) += delta(param_i); 
                }     
                if (constrain_parameter == 0){
                    current_parameters(param_j) += delta(param_j);
                }else if (constrain_parameter == 1 && param_j != 6){
                   current_parameters(param_j) += delta(param_j);
                }else if (constrain_parameter == 2 && param_j != 7){ 
                   current_parameters(param_j) += delta(param_j); 
                }                       
                    //current_parameters(i) += delta(i);
                   // current_parameters(param_i) += delta(param_i);
                   // current_parameters(param_i) += delta(param_j);     
                    int current_param_array[8] = { 0 };
                    current_param_array[i] = 1;
                    current_param_array[param_i] = 1;  
                    current_param_array[param_j] = 1;
                    double current_loglik = calculate_loglikelihood_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);                       
                    if(current_loglik > best_loglik){
                       best_loglik = current_loglik;
                             for(int j = 0 ; j < 8 ; j++){
                                param_array[j] = current_param_array[j];
                            }
                    }
               }
           } 
        }
           
         for(int i = 0; i < 8; i++){
            for(int param_a = i+1; param_a < 8 ; param_a++){
                for(int param_b = param_a + 1; param_b < 8 ; param_b++){
                    for(int param_c = param_b + 1; param_c < 8; param_c++){
                        Eigen::VectorXd current_parameters = parameters;
                        
                 if (constrain_parameter == 0){
                    current_parameters(i) += delta(i);
                }else if (constrain_parameter == 1 && i != 6){
                   current_parameters(i) += delta(i);
                }else if (constrain_parameter == 2 && i != 7){ 
                   current_parameters(i) += delta(i); 
                } 
                 if (constrain_parameter == 0){
                    current_parameters(param_a) += delta(param_a);
                }else if (constrain_parameter == 1 && param_a != 6){
                   current_parameters(param_a) += delta(param_a);
                }else if (constrain_parameter == 2 && param_a != 7){ 
                   current_parameters(param_a) += delta(param_a); 
                }     
                if (constrain_parameter == 0){
                    current_parameters(param_b) += delta(param_b);
                }else if (constrain_parameter == 1 && param_b != 6){
                   current_parameters(param_b) += delta(param_b);
                }else if (constrain_parameter == 2 && param_b != 7){ 
                   current_parameters(param_b) += delta(param_b); 
                } 
                if (constrain_parameter == 0){
                    current_parameters(param_c) += delta(param_c);
                }else if (constrain_parameter == 1 && param_c != 6){
                   current_parameters(param_c) += delta(param_c);
                }else if (constrain_parameter == 2 && param_c != 7){ 
                   current_parameters(param_c) += delta(param_c); 
                }                                                         
                      //  current_parameters(i) += delta(i);
                     //   current_parameters(param_a) += delta(param_a);
                      //  current_parameters(param_b) += delta(param_b);
                     //   current_parameters(param_c) += delta(param_c);
                        int current_param_array[8] = { 0 };
                        current_param_array[param_a] = 1;
                        current_param_array[param_b] = 1;
                        current_param_array[param_c] = 1;
                        current_param_array[i] = 1;
                        double current_loglik = calculate_loglikelihood_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);
                         if(current_loglik > best_loglik){
                             best_loglik = current_loglik;
                             for(int j = 0 ; j < 8 ; j++){
                                param_array[j] = current_param_array[j];
                            }
                         }                        
                        
                    }
               }
           }                                
       }
            
            
       
       
         for(int i = 0; i < 8; i++){
            for(int param_a = i+1; param_a < 8 ; param_a++){
                for(int param_b = param_a + 1; param_b < 8 ; param_b++){
                    for(int param_c = param_b + 1; param_c < 8; param_c++){
                        for(int param_d = param_c + 1; param_d < 8; param_d++){
                        Eigen::VectorXd current_parameters = parameters;
                 if (constrain_parameter == 0){
                    current_parameters(i) += delta(i);
                }else if (constrain_parameter == 1 && i != 6){
                   current_parameters(i) += delta(i);
                }else if (constrain_parameter == 2 && i != 7){ 
                   current_parameters(i) += delta(i); 
                } 
                 if (constrain_parameter == 0){
                    current_parameters(param_a) += delta(param_a);
                }else if (constrain_parameter == 1 && param_a != 6){
                   current_parameters(param_a) += delta(param_a);
                }else if (constrain_parameter == 2 && param_a != 7){ 
                   current_parameters(param_a) += delta(param_a); 
                }     
                if (constrain_parameter == 0){
                    current_parameters(param_b) += delta(param_b);
                }else if (constrain_parameter == 1 && param_b != 6){
                   current_parameters(param_b) += delta(param_b);
                }else if (constrain_parameter == 2 && param_b != 7){ 
                   current_parameters(param_b) += delta(param_b); 
                } 
                if (constrain_parameter == 0){
                    current_parameters(param_c) += delta(param_c);
                }else if (constrain_parameter == 1 && param_c != 6){
                   current_parameters(param_c) += delta(param_c);
                }else if (constrain_parameter == 2 && param_c != 7){ 
                   current_parameters(param_c) += delta(param_c); 
                }
                if (constrain_parameter == 0){
                    current_parameters(param_d) += delta(param_d);
                }else if (constrain_parameter == 1 && param_d != 6){
                   current_parameters(param_d) += delta(param_d);
                }else if (constrain_parameter == 2 && param_d != 7){ 
                   current_parameters(param_d) += delta(param_d); 
                }                                           
                       // current_parameters(i) += delta(i);
                      //  current_parameters(param_a) += delta(param_a);
                       // current_parameters(param_b) += delta(param_b);
                      //  current_parameters(param_c) += delta(param_c);
                     //   current_parameters(param_d) += delta(param_d);
                        int current_param_array[8] = { 0 };
                        current_param_array[param_a] = 1;
                        current_param_array[param_b] = 1;
                        current_param_array[param_c] = 1;
                        current_param_array[param_d] = 1;
                        current_param_array[i] = 1;
                        double current_loglik = calculate_loglikelihood_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);
                         if(current_loglik > best_loglik){
                             best_loglik = current_loglik;
                             for(int j = 0 ; j < 8 ; j++){
                                param_array[j] = current_param_array[j];
                            }
                         } 
                    }                       
               }
           }                                
         }
        }    
         for(int i = 0; i < 8; i++){
            for(int param_a = i+1; param_a < 8 ; param_a++){
                for(int param_b = param_a + 1; param_b < 8 ; param_b++){
                    for(int param_c = param_b + 1; param_c < 8; param_c++){
                        for(int param_d = param_c + 1; param_d < 8; param_d++){
                            for(int param_e = param_d + 1; param_e < 8; param_e++){
                        Eigen::VectorXd current_parameters = parameters;
                        
                 if (constrain_parameter == 0){
                    current_parameters(i) += delta(i);
                }else if (constrain_parameter == 1 && i != 6){
                   current_parameters(i) += delta(i);
                }else if (constrain_parameter == 2 && i != 7){ 
                   current_parameters(i) += delta(i); 
                } 
                 if (constrain_parameter == 0){
                    current_parameters(param_a) += delta(param_a);
                }else if (constrain_parameter == 1 && param_a != 6){
                   current_parameters(param_a) += delta(param_a);
                }else if (constrain_parameter == 2 && param_a != 7){ 
                   current_parameters(param_a) += delta(param_a); 
                }     
                if (constrain_parameter == 0){
                    current_parameters(param_b) += delta(param_b);
                }else if (constrain_parameter == 1 && param_b != 6){
                   current_parameters(param_b) += delta(param_b);
                }else if (constrain_parameter == 2 && param_b != 7){ 
                   current_parameters(param_b) += delta(param_b); 
                } 
                if (constrain_parameter == 0){
                    current_parameters(param_c) += delta(param_c);
                }else if (constrain_parameter == 1 && param_c != 6){
                   current_parameters(param_c) += delta(param_c);
                }else if (constrain_parameter == 2 && param_c != 7){ 
                   current_parameters(param_c) += delta(param_c); 
                }
                if (constrain_parameter == 0){
                    current_parameters(param_d) += delta(param_d);
                }else if (constrain_parameter == 1 && param_d != 6){
                   current_parameters(param_d) += delta(param_d);
                }else if (constrain_parameter == 2 && param_d != 7){ 
                   current_parameters(param_d) += delta(param_d); 
                }  
                if (constrain_parameter == 0){
                    current_parameters(param_e) += delta(param_e);
                }else if (constrain_parameter == 1 && param_e != 6){
                   current_parameters(param_e) += delta(param_e);
                }else if (constrain_parameter == 2 && param_e != 7){ 
                   current_parameters(param_e) += delta(param_e); 
                }                                         
                     //   current_parameters(i) += delta(i);
                    //    current_parameters(param_a) += delta(param_a);
                   //     current_parameters(param_b) += delta(param_b);
                    //    current_parameters(param_c) += delta(param_c);
                   //     current_parameters(param_d) += delta(param_d);
                 //       current_parameters(param_e) += delta(param_e);
                        int current_param_array[8] = { 0 };
                        current_param_array[param_a] = 1;
                        current_param_array[param_b] = 1;
                        current_param_array[param_c] = 1;
                        current_param_array[param_d] = 1;
                        current_param_array[param_e] = 1;
                        current_param_array[i] = 1;
                        double current_loglik = calculate_loglikelihood_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);
                         if(current_loglik > best_loglik){
                             best_loglik = current_loglik;
                             for(int j = 0 ; j < 8 ; j++){
                                param_array[j] = current_param_array[j];
                            }
                         } 
                    }                       
               }
           }                                
         }
        }             
       }       
       
         for(int i = 0; i < 8; i++){
            for(int param_a = i+1; param_a < 8 ; param_a++){
                for(int param_b = param_a + 1; param_b < 8 ; param_b++){
                    for(int param_c = param_b + 1; param_c < 8; param_c++){
                        for(int param_d = param_c + 1; param_d < 8; param_d++){
                            for(int param_e = param_d + 1; param_e < 8; param_e++){
                                 for(int param_f = param_e + 1; param_f < 8; param_f++){
                        Eigen::VectorXd current_parameters = parameters;
                 if (constrain_parameter == 0){
                    current_parameters(i) += delta(i);
                }else if (constrain_parameter == 1 && i != 6){
                   current_parameters(i) += delta(i);
                }else if (constrain_parameter == 2 && i != 7){ 
                   current_parameters(i) += delta(i); 
                } 
                 if (constrain_parameter == 0){
                    current_parameters(param_a) += delta(param_a);
                }else if (constrain_parameter == 1 && param_a != 6){
                   current_parameters(param_a) += delta(param_a);
                }else if (constrain_parameter == 2 && param_a != 7){ 
                   current_parameters(param_a) += delta(param_a); 
                }     
                if (constrain_parameter == 0){
                    current_parameters(param_b) += delta(param_b);
                }else if (constrain_parameter == 1 && param_b != 6){
                   current_parameters(param_b) += delta(param_b);
                }else if (constrain_parameter == 2 && param_b != 7){ 
                   current_parameters(param_b) += delta(param_b); 
                } 
                if (constrain_parameter == 0){
                    current_parameters(param_c) += delta(param_c);
                }else if (constrain_parameter == 1 && param_c != 6){
                   current_parameters(param_c) += delta(param_c);
                }else if (constrain_parameter == 2 && param_c != 7){ 
                   current_parameters(param_c) += delta(param_c); 
                }
                if (constrain_parameter == 0){
                    current_parameters(param_d) += delta(param_d);
                }else if (constrain_parameter == 1 && param_d != 6){
                   current_parameters(param_d) += delta(param_d);
                }else if (constrain_parameter == 2 && param_d != 7){ 
                   current_parameters(param_d) += delta(param_d); 
                }  
                if (constrain_parameter == 0){
                    current_parameters(param_e) += delta(param_e);
                }else if (constrain_parameter == 1 && param_e != 6){
                   current_parameters(param_e) += delta(param_e);
                }else if (constrain_parameter == 2 && param_e != 7){ 
                   current_parameters(param_e) += delta(param_e); 
                }  
                if (constrain_parameter == 0){
                    current_parameters(param_f) += delta(param_f);
                }else if (constrain_parameter == 1 && param_f != 6){
                   current_parameters(param_f) += delta(param_f);
                }else if (constrain_parameter == 2 && param_f != 7){ 
                   current_parameters(param_f) += delta(param_f); 
                }                                         
                        
                     //   current_parameters(i) += delta(i);
                    //    current_parameters(param_a) += delta(param_a);
                    //    current_parameters(param_b) += delta(param_b);
                    //    current_parameters(param_c) += delta(param_c);
                    //    current_parameters(param_d) += delta(param_d);
                    //    current_parameters(param_e) += delta(param_e);
                   //     current_parameters(param_f) += delta(param_f);
                        int current_param_array[8] = { 0 };
                        current_param_array[param_a] = 1;
                        current_param_array[param_b] = 1;
                        current_param_array[param_c] = 1;
                        current_param_array[param_d] = 1;
                        current_param_array[param_e] = 1;
                        current_param_array[param_f] = 1;
                        current_param_array[i] = 1;
                        double current_loglik = calculate_loglikelihood_param(trait_one, trait_two, Ones, eigenvalues, current_parameters);
                         if(current_loglik > best_loglik){
                             best_loglik = current_loglik;
                             for(int j = 0 ; j < 8 ; j++){
                                param_array[j] = current_param_array[j];
                            }
                           
                         } 
                    }                      
               }
           }                                
         }
        }             
       }
      }*/

              for(int i = 0 ; i < 8; i++){
                if (param_array[i] == 1){
                    if(constrain_parameter == 1 && i == 6)
                        continue;
                    if(constrain_parameter == 2 && i == 7)
                        continue;
                double relative_change = abs(delta(i)/parameters(i)); 
                if(relative_change != relative_change) 
                    continue;            
                parameters(i) += delta(i);
                if ((max_delta < relative_change) || (max_delta == -1 && relative_change != 0))
                       max_delta = relative_change;
            }
                
      }  
      
     
     double h2r_one = calculate_constraint(parameters(2)); 
     double h2r_two = calculate_constraint(parameters(3)); 
     double rhog = calculate_rho(parameters(6)); 
     double rhoe = calculate_rho(parameters(7));
   /*  double test_h2r_one = h2r_one;
     double test_t_one = parameters(2); 
     double test_h2r_two = h2r_one;
     double test_t_two = parameters(3);
     double test_rhog = rhog;
     double test_i_one = parameters(6);   
     double test_rhoe = rhoe;
     double test_i_two = parameters(7); 
     bool run_new_loglik_test = false;                              
        if(h2r_one >= 0.9999){
                test_h2r_one = 1.0;
                test_t_one = 1000;
                run_new_loglik_test = true;
               
         }
         
         if(h2r_one <= 1e-05){
                test_h2r_one = 0.0;
                test_t_one  = 0.001;
                run_new_loglik_test = true;
         }
         
         if(h2r_two >= 0.9999){
                test_h2r_two = 1.0;
                test_t_two = 1000;
                run_new_loglik_test = true;
         }
                           
         if(h2r_two <= 1e-05){
                test_h2r_two = 0.0;
               test_t_two = 0.001;
               run_new_loglik_test = true;
         } 
                  
        if (rhog >= 0.99 && constrain_parameter != 1){
              test_rhog = 1.0;
              test_i_one = 1000000000;
              run_new_loglik_test = true;
        }  
        
        if (rhog <= -0.99 && constrain_parameter != 1){
              rhog = -1.0;
              test_i_one  = -1000000000;
              run_new_loglik_test = true;
        }  
        
        if (rhoe >= 0.99 && constrain_parameter != 2){
              test_rhoe = 1.0;
              test_i_two = 1000000000;
              run_new_loglik_test = true;
        }  
        
        if (rhoe <= -0.99 && constrain_parameter != 2){
              test_rhoe = -1.0;      
              test_i_two  = -1000000000;
              run_new_loglik_test = true;
        }
        
        if (rhog >= 0.999999){
            test_rhog = 1.0;
            test_i_one = 1000000000000;
            run_new_loglik_test = true;
        }
        if (rhog <= -0.999999){
            test_rhog = -1.0;
            test_i_one = -1000000000000;
            run_new_loglik_test = true;
        } 
        if (rhoe >= 0.999999){
            test_rhoe = 1.0;
            test_i_two = 1000000000000;
            run_new_loglik_test = true;
        }
        if (rhoe <= -0.999999){
            test_rhoe = -1.0;
            test_i_two = -1000000000000;
            run_new_loglik_test = true;
        }      */         
         last_loglik = loglik; 
         loglik =  calculate_loglikelihood_param(trait_one, trait_two, Ones, eigenvalues, parameters);                                                            

        /* double test_loglik = loglik;
         if (run_new_loglik_test){
            test_loglik =  calculate_loglikelihood( trait_one, trait_two, Ones, eigenvalues, parameters(0) ,  parameters(1) , \
                                              parameters(4) ,parameters(5) , h2r_one,  h2r_two, test_rhog, test_rhoe);

        }
                                                                                                                                            
       if(test_loglik >= loglik && run_new_loglik_test && test_loglik == test_loglik && loglik == loglik){
            
            if(test_h2r_one != h2r_one){
                h2r_one = test_h2r_one;
                parameters(2) = test_t_one;
            }
            if(test_h2r_two != h2r_two){
                h2r_two = test_h2r_two;
                parameters(3) = test_t_two;
            }*
           
             rhog = test_rhog;
             parameters(6) = test_i_one;
            
          
            parameters(7) = test_i_two;
           rhoe = test_rhoe;
         
           loglik = test_loglik;                        
       } */                                                                                                                                                                                                                                                                                                          
      if (debug == true){

        cout << "mean one: " << parameters(0) << " mean two: " << parameters(1) << endl;
        cout << "sd one: " << parameters(4) << " sd two: " << parameters(5) << endl;
        cout << "t one: " << parameters(2) << " h2r one: " << h2r_one << endl;
        cout << "t two: " << parameters(3) << " h2r two: " << h2r_two << endl;
        cout << "i one: " << parameters(6) << " rhog: " << rhog << endl;
        cout << "i two: " << parameters(7) << " rhoe: " << rhoe << endl;
        cout << "loglik: " << loglik << endl;

     } 

        iteration++; 
    }while(loglik == loglik && delta == delta && iteration < MAX_ITERATIONS  && (abs((last_loglik - loglik)/last_loglik) >= MAX_LOGLIK_ERROR || max_delta >= MAX_DELTA_ERROR));
    #pragma omp critical 
    {
    if(loglik == loglik && delta == delta && (loglik > best_loglikelihood || first_loglik == false)){
        best_h = h;
        best_loglikelihood = loglik;
        final_parameters = parameters;
        first_loglik = true;
        iteration_count = iteration;
     }
    }
  }
  if(first_loglik == false){
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
/*
    Eigen::VectorXd residual_one = trait_one - covariate_matrix*(covariate_matrix.transpose()*covariate_matrix).inverse()*covariate_matrix.transpose()*trait_one;
    Eigen::VectorXd residual_two = trait_two - covariate_matrix*(covariate_matrix.transpose()*covariate_matrix).inverse()*covariate_matrix.transpose()*trait_two;
    Eigen::MatrixXd aux = Eigen::MatrixXd::Ones(residual_one.rows(), 2);
    aux.col(1) = eigenvalues;
    Eigen::VectorXd theta_one = (aux.transpose()*aux).inverse()*aux.transpose()*(residual_one.cwiseAbs2());
    Eigen::VectorXd theta_two = (aux.transpose()*aux).inverse()*aux.transpose()*(residual_two.cwiseAbs2());
    double var_one = theta_one(0) + theta_one(1);
    double var_two = theta_two(0) + theta_two(1);
    double h2r_one = theta_one(1)/var_one;
    double h2r_two = theta_two(1)/var_two;
    Eigen::MatrixXd aux_two = aux;
    aux_two.col(0) = aux_two.col(0)*sqrt((1.0-h2r_one)*(1.0-h2r_two)*var_one*var_two);
    aux_two.col(1) = aux_two.col(1)*sqrt((h2r_one)*(h2r_two)*var_one*var_two);
    Eigen::VectorXd theta_three = (aux_two.transpose()*aux_two).inverse()*aux_two.transpose()*(residual_one.cwiseProduct(residual_two));
    if(theta_three(0) > 1.0)
        theta_three(0) = 1.0;
    if(theta_three(1) > 1.0)
        theta_three(1) = 1.0;
    if(theta_three(0) < -1.0)
        theta_three(0) = -1.0;
    if(theta_three(1) < -1.0)
        theta_three(1) = -1.0;
    cout << theta_three << endl << endl;
    cout << theta_one(1)/(theta_one(1) + theta_one(0)) << endl << endl;
    cout << theta_two(1)/(theta_two(1) + theta_two(0)) << endl << endl;
    double rhop_test = theta_three(0)*sqrt((1.0-h2r_one)*(1.0-h2r_two)) +  theta_three(1)*sqrt((h2r_one)*(h2r_two));                         
    Eigen::MatrixXd omega_one = (aux*theta_one).cwiseInverse().cwiseAbs2().asDiagonal()/(1.0-rhop*rhop);
    Eigen::MatrixXd omega_two = (aux*theta_two).cwiseInverse().cwiseAbs2().asDiagonal();
    
    theta_one = (aux.transpose()*omega_one*aux).inverse()*aux.transpose()*omega_one*(residual_one.cwiseAbs2());
    theta_two = (aux.transpose()*omega_two*aux).inverse()*aux.transpose()*omega_two*(residual_two.cwiseAbs2());
    var_one = theta_one(0) + theta_one(1);
    var_two = theta_two(0) + theta_two(1);
    h2r_one = theta_one(1)/var_one;
    h2r_two = theta_two(1)/var_two;
    aux_two = aux;
    aux_two.col(0) = aux_two.col(0)*sqrt((1.0-h2r_one)*(1.0-h2r_two)*var_one*var_two);
    aux_two.col(1) = aux_two.col(1)*sqrt((h2r_one)*(h2r_two)*var_one*var_two);
    Eigen::MatrixXd omega_three = (aux_two*theta_three).cwiseInverse().cwiseAbs2().asDiagonal();    
    theta_three = (aux_two.transpose()*omega_three*aux_two).inverse()*aux_two.transpose()*omega_three*(residual_one.cwiseProduct(residual_two));

    double rhog_test = theta_three(1); 
    double rhoe_test = theta_three(0);   
    cout << "h2r one: " << h2r_one << endl;
    cout << "h2r two: " << h2r_two  << endl;
    cout << "rhog: " << rhog_test << endl;
    cout << "rhoe: " << rhoe_test << endl;
    cout << "rhop: " << rhog_test*sqrt(h2r_one*h2r_two) + rhoe_test*sqrt((1.0-h2r_one)*(1.0-h2r_two)) << endl;

    return 0;*/

    Eigen::VectorXd final_parameters(8);
    Eigen::VectorXd final_beta(covariate_matrix.cols()*2);
   // double main_best_h;
    int main_iteration_count;
    Eigen::VectorXd standard_errors;
    double main_loglik = calculate_bivariate_model_three(trait_one,  trait_two,  covariate_matrix, eigenvalues,  final_parameters, final_beta, standard_errors, main_best_h, main_iteration_count, 0,  debug); ;// calculate_bivariate_model_two(trait_one,  trait_two,  S_ones, eigenvalues,  final_parameters,  main_best_h, main_iteration_count, 0,  debug); 
    if(main_loglik != main_loglik ){
        Solar_Eval(interp, "loglike set 0");
        const char * error = "Convergence of this model could not be achieved";
        delete file_data;
        return error;
    }
    string loglik_command = "loglike set " + to_string(main_loglik);
    Solar_Eval(interp, loglik_command.c_str());   
    if(display_pvalues){
        Eigen::VectorXd rhog_parameters(8);
        Eigen::VectorXd rhog_beta(covariate_matrix.cols()*2);
        int iteration_count;
        double rhog_best_h = main_best_h;
        Eigen::VectorXd rhog_standard_errors;
        double rhog_loglik  = calculate_bivariate_model_three(trait_one,  trait_two,  covariate_matrix, eigenvalues,  rhog_parameters, rhog_beta, rhog_standard_errors, rhog_best_h, main_iteration_count, 1,  debug); ;//calculate_bivariate_model_two(trait_one,  trait_two,  S_ones, eigenvalues,  rhog_parameters,  rhog_best_h, iteration_count, 1,  debug);
        
        
        
        Eigen::VectorXd rhoe_parameters(8);
        Eigen::VectorXd rhoe_beta(covariate_matrix.cols()*2);
        double rhoe_best_h = main_best_h;
        Eigen::VectorXd rhoe_standard_errors;
        double rhoe_loglik  = calculate_bivariate_model_three(trait_one,  trait_two,  covariate_matrix, eigenvalues,  rhoe_parameters, rhoe_beta, rhoe_standard_errors, rhoe_best_h, main_iteration_count, 2,  debug); ;//calculate_bivariate_model_two(trait_one,  trait_two,  S_ones, eigenvalues,  rhoe_parameters,  rhoe_best_h, iteration_count, 2,  debug);

       
        double rhog_pvalue = 1.0;
        if(rhog_loglik == rhog_loglik){
            rhog_pvalue = 2.0*chicdf(2.0*(main_loglik - rhog_loglik), 1); 
        }
        double rhoe_pvalue = 1.0;
        if(rhoe_loglik == rhoe_loglik){
            rhoe_pvalue = 2.0*chicdf(2.0*(main_loglik - rhoe_loglik), 1);
        }
        const double rhop = calculate_rho(final_parameters(6))*sqrt(calculate_constraint(final_parameters(2)))*\
             sqrt(calculate_constraint(final_parameters(3))) + calculate_rho(final_parameters(7))*\
             sqrt(1.0-calculate_constraint(final_parameters(2)))*sqrt(1.0-calculate_constraint(final_parameters(3)));
        final_parameters(2) = calculate_constraint(final_parameters(2));
        final_parameters(3) = calculate_constraint(final_parameters(3)); 
        final_parameters(6) = calculate_rho(final_parameters(6));
        final_parameters(7) = calculate_rho(final_parameters(7));               
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
        
        cout << setw(largest_name + 4) << parameter_names[6 +covariate_matrix.cols() - 1] << setw(20) << final_parameters(0) << setw(20) << standard_errors(6 + covariate_matrix.cols() - 1) << endl;
        cout << setw(largest_name + 4) << parameter_names[6 + covariate_matrix.cols()*2 - 1] << setw(20) << final_parameters(1) << setw(20) << standard_errors(6 + 2*covariate_matrix.cols() - 1) << endl;
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
        output_stream << parameter_names[6 +covariate_matrix.cols() - 1] << "," << final_parameters(0) << "," << standard_errors(6 + covariate_matrix.cols() - 1) << endl;
        output_stream << parameter_names[6 + covariate_matrix.cols()*2  - 1] << "," << final_parameters(1) << "," << standard_errors(6 + 2*covariate_matrix.cols() - 1) << endl;
        for(int i = 0 ; i < covariate_matrix.cols() - 1; i++){
            output_stream << parameter_names[6 + i] << "," << final_beta(i) << ","<< standard_errors(6 + i) << endl;
            output_stream << parameter_names[6 + covariate_matrix.cols() + i] << "," << final_beta(i + covariate_matrix.cols()) << ","<< standard_errors(6 + covariate_matrix.cols() + i) << endl;
        }

        for(int i = 0 ; i < 6; i++){
            output_stream << parameter_names[i] << "," << final_parameters(2 + i)<< ","<< standard_errors(i) << endl;
             
        } 
        output_stream << "rhop," << rhop << ",\n";
        output_stream << "loglik," << main_loglik << ",\n";
        output_stream << "rhog p-value," << rhog_pvalue << ",\n";
        output_stream << "rhog loglik," << rhog_loglik << ",\n";
        output_stream << "rhoe p-value," << rhoe_pvalue << ",\n"; 
        output_stream << "rhoe loglik," << rhoe_loglik << ",\n"; 
        output_stream.close();
        /*output_stream << parameter_names[6 + covariate_matrix.cols() - 1] <<  "," << final_parameters(0) << standard_errors(6 + covariate_matrix.cols() - 1)endl;
        output_stream << trait_list[1] <<  "-mean," << final_parameters(1) << endl;
        for(int i = 0 ; i < covariate_matrix.cols() -1 ; i++){
            output_stream << trait_list[0] <<  "-b" << cov_list[i] << ": " << final_beta(i) << endl;
            output_stream << trait_list[1] <<  "-b" << cov_list[i] << ": "  << final_beta(covariate_matrix.cols() + i) << endl;
        }        
        output_stream << trait_list[0] <<  "-sd," << abs(final_parameters(4)) << endl;
        output_stream << trait_list[1] <<  "-sd," << abs(final_parameters(5)) << endl;
        output_stream << trait_list[0] <<  "-h2r," << calculate_constraint(final_parameters(2)) << endl;
        output_stream << trait_list[1] <<  "-h2r," << calculate_constraint(final_parameters(3)) << endl;
        output_stream << "rhog," << calculate_rho(final_parameters(6)) << endl;
        output_stream << "rhoe," << calculate_rho(final_parameters(7)) << endl;
        output_stream << "rhop," << rhop << endl;
        output_stream << "rhog-p-value," << rhog_pvalue<< endl;
        output_stream << "rhoe-p-value," << rhoe_pvalue<< endl;
        output_stream << "loglik," << main_loglik << endl;
        output_stream.close();*/
       // Eigen::VectorXd errors = calculate_errors(trait_one, trait_two, S_ones,eigenvalues, final_parameters(0), final_parameters(1),\
                                        final_parameters(2), final_parameters(3), final_parameters(4), final_parameters(5), final_parameters(6), final_parameters(7), main_best_h);
        //cout << endl;   
        //Solar_Eval(interp, solar_command.c_str());                              
        /*Eigen::VectorXd gradient(8);
        calculate_gradient(gradient, trait_one, trait_two,S_ones,  eigenvalues, final_parameters(0), final_parameters(1), final_parameters(2),final_parameters(3),\
                            final_parameters(4), final_parameters(5), final_parameters(6), final_parameters(7), main_best_h); 
        cout << "Gradient of Parameters\n";
        cout <<  trait_list[0] <<  "-mean: " << gradient(0) << endl;
         cout <<  trait_list[1] <<  "-mean: " << gradient(1) << endl;
         cout <<  trait_list[0] <<  "-h2r: " << gradient(2) << endl;
          cout <<  trait_list[1] <<  "-h2r: " << gradient(3) << endl;
          cout <<  trait_list[0] <<  "-sd: " << gradient(4) << endl;
        cout <<  trait_list[1] <<  "-sd: " << gradient(5) << endl;
        cout <<  "rhog: " << gradient(6) << endl;
        cout <<  "rhoe: " << gradient(7)<< endl;  */      
        
    }else{
        const double rhop = calculate_rho(final_parameters(6))*sqrt(calculate_constraint(final_parameters(2)))*\
             sqrt(calculate_constraint(final_parameters(3))) + calculate_rho(final_parameters(7))*\
             sqrt(1.0-calculate_constraint(final_parameters(2)))*sqrt(1.0-calculate_constraint(final_parameters(3)));
        final_parameters(2) = calculate_constraint(final_parameters(2));
        final_parameters(3) = calculate_constraint(final_parameters(3)); 
        final_parameters(6) = calculate_rho(final_parameters(6));
        final_parameters(7) = calculate_rho(final_parameters(7));  
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
        
        cout << setw(largest_name + 4) << parameter_names[6 +covariate_matrix.cols() - 1] << setw(20) << final_parameters(0) << setw(20) << standard_errors(6 + covariate_matrix.cols() - 1) << endl;
        cout << setw(largest_name + 4) << parameter_names[6 + 2*covariate_matrix.cols()  - 1] << setw(20) << final_parameters(1) << setw(20) << standard_errors(6 + 2*covariate_matrix.cols() - 1) << endl;
        for(int i = 0 ; i < covariate_matrix.cols() - 1; i++){
            cout << setw(largest_name + 4) << parameter_names[6 + i] << setw(20) << final_beta(i) << setw(20) << standard_errors(6 + i) << endl;
            cout << setw(largest_name + 4) << parameter_names[6 + covariate_matrix.cols() + i] << setw(20) << final_beta(i + covariate_matrix.cols()) << setw(20) << standard_errors(6 + covariate_matrix.cols() + i) << endl;
        }

        for(int i = 0 ; i < 6; i++){
            cout << setw(largest_name + 4) << parameter_names[i] << setw(20) << final_parameters(2 + i)<< setw(20) << standard_errors(i) << endl;
             
        }
        cout << setw(largest_name + 4) << "rhop" << setw(20) << rhop << endl;
        cout << setw(largest_name + 4) << "loglik" << setw(20) << main_loglik << endl;
 


        string output_filename = "gen_corr-" + trait_list[0] + "-" + trait_list[1] + ".out";
        ofstream output_stream(output_filename);
        output_stream << "Parameter,Value,Standard Error\n";
        output_stream << parameter_names[6 +covariate_matrix.cols() - 1] << "," << final_parameters(0) << "," << standard_errors(6 + covariate_matrix.cols() - 1) << endl;
        output_stream << parameter_names[6 + covariate_matrix.cols()*2 - 1] << "," << final_parameters(1) << "," << standard_errors(6 + 2*covariate_matrix.cols() - 1) << endl;
        for(int i = 0 ; i < covariate_matrix.cols() - 1; i++){
            output_stream << parameter_names[6 + i] << "," << final_beta(i) << "," << standard_errors(6 + i) << endl;
            output_stream << parameter_names[6 + covariate_matrix.cols() + i] << "," << final_beta(i + covariate_matrix.cols()) << ","<< standard_errors(6 + covariate_matrix.cols() + i) << endl;
        }

        for(int i = 0 ; i < 6; i++){
            output_stream << parameter_names[i] << "," << final_parameters(2 + i)<< ","<< standard_errors(i) << endl;
             
        } 
        output_stream << "rhop," << rhop << ",\n";
        output_stream << "loglik," << main_loglik << ",\n";
        output_stream.close();             
        //string solar_command = "set RHOP " + to_string(rhop);
       /* cout << "*******************************" << endl;
        cout << "*     Genetic Correlation     *" << endl;
        cout << "*******************************" << endl;                
        
        cout << "Pedigree: " << currentPed->filename() << endl;
        cout << "Phenotype: " << phenotype_filename << endl;
        cout << "Number of Subjects: " << eigenvalues.rows() << endl; 
        cout << "Total iterations: " << main_iteration_count << endl;
        cout << "Numerical Differentiation delta: " << main_best_h << endl << endl;
        cout << "Final parameters\n";
        cout << trait_list[0] <<  "-mean: " << final_parameters(0) << " " << trait_list[1] << "-mean: "  << final_parameters(1) << endl;
        for(int i = 0 ; i < covariate_matrix.cols() -1 ; i++){
            cout << trait_list[0] <<  "-b" << cov_list[i] << ": " << final_beta(i) << " " << trait_list[1] <<  "-b" << cov_list[i] << ": "  << final_beta(covariate_matrix.cols() + i) << endl;
        }        
        cout << trait_list[0] << "-h2r: " << calculate_constraint(final_parameters(2)) << " " << trait_list[1] << "-h2r: " << calculate_constraint(final_parameters(3)) << endl;
        cout << trait_list[0] << "-sd: "  << abs(final_parameters(4)) << " " << trait_list[1] << "-sd: " << abs(final_parameters(5)) << endl;
        cout << "rhog: " << calculate_rho(final_parameters(6)) << " rhoe: " << calculate_rho(final_parameters(7))\
             << " rhop: " << rhop << endl; 
        cout << "loglik: " << main_loglik << endl; 
        string output_filename = "gen_corr-" + trait_list[0] + "-" + trait_list[1] + ".out";
        ofstream output_stream(output_filename);
        output_stream << trait_list[0] <<  "-mean," << final_parameters(0) << endl;
        output_stream << trait_list[1] <<  "-mean," << final_parameters(1) << endl;
        for(int i = 0 ; i < covariate_matrix.cols() -1 ; i++){
            output_stream << trait_list[0] <<  "-b" << cov_list[i] << ": " << final_beta(i) << endl;
            output_stream << trait_list[1] <<  "-b" << cov_list[i] << ": "  << final_beta(covariate_matrix.cols() + i) << endl;
        }        
        output_stream << trait_list[0] <<  "-sd," << abs(final_parameters(4)) << endl;
        output_stream << trait_list[1] <<  "-sd," << abs(final_parameters(5)) << endl;
        output_stream << trait_list[0] <<  "-h2r," << calculate_constraint(final_parameters(2)) << endl;
        output_stream << trait_list[1] <<  "-h2r," << calculate_constraint(final_parameters(3)) << endl;
        output_stream << "rhog," << calculate_rho(final_parameters(6)) << endl;
        output_stream << "rhoe," << calculate_rho(final_parameters(7)) << endl;
        output_stream << "rhop," << rhop << endl;
        output_stream << "loglik," << main_loglik << endl;
        output_stream.close();*/
        //cout << endl;
      //  Solar_Eval(interp, solar_command.c_str()); 
       /* Eigen::VectorXd gradient(8);
        calculate_gradient(gradient, trait_one, trait_two,S_ones,  eigenvalues, final_parameters(0), final_parameters(1), final_parameters(2),final_parameters(3),\
                            final_parameters(4), final_parameters(5), final_parameters(6), final_parameters(7), main_best_h) ;
        cout << "Gradient of Parameters\n";
        cout <<  trait_list[0] <<  "-mean: " << gradient(0) << endl;
         cout <<  trait_list[1] <<  "-mean: " << gradient(1) << endl;
         cout <<  trait_list[0] <<  "-h2r: " << gradient(2) << endl;
          cout <<  trait_list[1] <<  "-h2r: " << gradient(3) << endl;
          cout <<  trait_list[0] <<  "-sd: " << gradient(4) << endl;
        cout <<  trait_list[1] <<  "-sd: " << gradient(5) << endl;
        cout <<  "rhog: " << gradient(6) << endl;
        cout <<  "rhoe: " << gradient(7)<< endl; 
        cout << endl;
        Eigen::VectorXd errors = calculate_errors(trait_one, trait_two, S_ones, eigenvalues, final_parameters(0), final_parameters(1), final_parameters(2),final_parameters(3),\
                            final_parameters(4), final_parameters(5), final_parameters(6), final_parameters(7), 0.1) ;
        cout << "Standard Errors\n";
          cout << errors << endl;*/
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
            
            
static double calculate_bivariate_model_four(Eigen::VectorXd trait_one, Eigen::VectorXd trait_two, Eigen::MatrixXd covariate_matrix, Eigen::VectorXd eigenvalues, Eigen::VectorXd & final_parameters, Eigen::VectorXd & final_beta, double  & best_h, int & iteration_count, int constrain_parameter = 0, bool debug = false){


    double best_loglikelihood;
    bool first_loglik = false;
    double h;
    omp_set_num_threads(10);
    const int index_map[6] = {2, 3, 4, 5, 6, 7};
    Eigen::VectorXd beta(covariate_matrix.cols()*2);
    Eigen::VectorXd parameters(8);
    parameters(2) = 0.35;
    parameters(3) = 0.35;
    parameters(4) = sqrt(((trait_one.array() - trait_one.mean()).matrix().squaredNorm()/trait_one.rows()));
    parameters(5) = sqrt(((trait_two.array() - trait_two.mean())).matrix().squaredNorm()/trait_two.rows());
    parameters(6) = 0.025;
    parameters(7) = 0.025;
    Eigen::VectorXd delta(6);
    Eigen::MatrixXd hessian(6, 6);
    Eigen::VectorXd gradient(6);    
    calculate_beta(trait_one, trait_two, covariate_matrix, eigenvalues, parameters ,beta);
    genetic_correlation_calculate_hessian_and_gradient_multithread_four(hessian,  gradient, trait_one, trait_two,\
                                                               covariate_matrix,  eigenvalues,  parameters, beta);
    Eigen::VectorXd test_vector = (-hessian.inverse()*gradient.cwiseAbs2()).cwiseAbs();
    delta = -hessian.inverse()*gradient;
    double largest_test = 0.0;
    for(int i = 0 ; i < 6; i++){
        if(test_vector(i) > largest_test){
            largest_test = test_vector(i);
        }
    }
    int iteration = 0;

    double last_loglik = 0.0;
    double loglik = 0.0;
    double max_delta;
    while(largest_test > STOP_CRITERIA && delta == delta && iteration <MAX_ITERATIONS && loglik == loglik){


        delta = -hessian.inverse()*gradient;
      //  cout << "delta: " << delta << endl;                
        //max_delta = -1;

        int param_array[6] = {1, 1, 1, 1, 1, 1};

        if(constrain_parameter == 1)
            delta(4) = 0.0;
        if(constrain_parameter == 2)
            delta(5) = 0.0;

        Eigen::VectorXd test_parameters = parameters;

        for(int i = 0 ; i < 6; i++)
            test_parameters(index_map[i]) += T*delta(i);

        Eigen::VectorXd test_beta = beta;
        calculate_beta(trait_one, trait_two, covariate_matrix, eigenvalues, test_parameters ,test_beta);
        double best_loglik = calculate_loglikelihood_param_two(trait_one, trait_two, covariate_matrix, eigenvalues, test_parameters, test_beta);
        Eigen::VectorXd best_parameters = test_parameters;
        Eigen::VectorXd best_beta = test_beta;
        for(int param_index = 1; param_index < 63; param_index++){
            int current_param_array[6] = {0};
            Eigen::VectorXd current_parameters = parameters;
            Eigen::VectorXd current_beta = beta;
            if((param_index & 1) == 0){
                current_parameters(2) += T*delta(0);
                current_param_array[0] = 1;
            }
            if((param_index & 2) == 0){
                current_parameters(3) += T*delta(1);
                current_param_array[1] = 1;
            }
            if((param_index & 4) == 0){
                current_parameters(4) += T*delta(2);
                current_param_array[2] = 1;
            }
            if((param_index & 8) == 0){
                current_parameters(5) += T*delta(3);
                current_param_array[3] = 1;
            }
            if((param_index & 16) == 0){
                current_parameters(6) += T*delta(4);
                current_param_array[4] = 1;
            }
            if((param_index & 32) == 0){
                current_parameters(7) += T*delta(5);
                current_param_array[5] = 1;
            }                        
            calculate_beta(trait_one, trait_two, covariate_matrix, eigenvalues, current_parameters, current_beta);
            double current_loglik = calculate_loglikelihood_param_two(trait_one, trait_two, covariate_matrix, eigenvalues, current_parameters, current_beta);
            if(current_loglik > best_loglik || best_loglik == 0.0){
                best_loglik = current_loglik;
                best_parameters = current_parameters;
                best_beta = current_beta;
                for(int j = 0; j < 6; j++){
                    param_array[j] = current_param_array[j];
                }
            }            
        }
        parameters = best_parameters;
        beta = best_beta;
        loglik = best_loglik;
/*
        for(int i = 0 ; i < 6; i++){
            if (param_array[i] == 1){
                if(constrain_parameter == 1 && i == 4)
                    continue;
                if(constrain_parameter == 2 && i == 5)
                    continue;
                parameters(index_map[i]) += delta(i);
     
            }
                
      }
        calculate_beta(trait_one, trait_two, covariate_matrix, eigenvalues, parameters, beta);
        loglik =  calculate_loglikelihood_param_two(trait_one, trait_two, covariate_matrix, eigenvalues, parameters, beta); */
        if(debug){
            cout << "Iteration: " << iteration << endl;
            cout << "h2 one: " << calculate_constraint(parameters(2)) << " " << " h2 two: " << calculate_constraint(parameters(3)) << endl;
            cout << "sd one: " << abs(parameters(4)) << " " << " sd two: " << abs(parameters(5)) << endl;
            cout << "rhog: " << calculate_rho(parameters(6)) << " " << " rhoe: " << calculate_rho(parameters(7)) << endl;
            cout << "loglik: " << loglik << endl;

        }
          
        genetic_correlation_calculate_hessian_and_gradient_multithread_four(hessian,  gradient, trait_one, trait_two,\
                                                               covariate_matrix,  eigenvalues,  parameters, beta);
        test_vector = (-hessian.inverse()*gradient.cwiseAbs2()).cwiseAbs();
        for(int i = 0 ; i < 6; i++){
            if(test_vector(i) > largest_test){
                largest_test = test_vector(i);
            }
        }    
        iteration++;                                                              
    }
    if(loglik == loglik && delta == delta && loglik == loglik && iteration <= MAX_ITERATIONS){
        best_loglikelihood = loglik;
        final_parameters = parameters;
        final_beta = beta;

        iteration_count = iteration;
        return best_loglikelihood;
    }else{
        return nan("");
    }
  
} 