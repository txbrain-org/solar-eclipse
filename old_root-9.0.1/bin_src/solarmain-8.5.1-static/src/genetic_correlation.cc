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
#define MAX_ITERATIONS 300
#define MAX_DELTA_ERROR 1e-07
#define MAX_LOGLIK_ERROR 1e-10
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
static double calculate_bivariate_model(Eigen::VectorXd trait_one, Eigen::VectorXd trait_two, Eigen::VectorXd Ones, Eigen::VectorXd eigenvalues, Eigen::VectorXd & final_parameters, double  & best_h, int & iteration_count, int constrain_parameter = 0, bool debug = false){

    double best_loglikelihood;
    bool first_loglik = false;
    double h;
    omp_set_num_threads(omp_get_max_threads());
   
#pragma omp parallel for private(h) shared(first_loglik)
    for(int i = 1; i <= 10 ; i++){
        if (i == 1){ 
         cout << "Using " << omp_get_num_threads() << " threads for computation\n"; 
        }
        h = 5*pow(10, -i);
        
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
static const char * calculate_genetic_correlation(Tcl_Interp * interp, vector<string> trait_list, const char* phenotype_filename, bool display_pvalues, bool debug, const char * evd_data_filename = 0){
    solar_mle_setup * file_data;
    try{
        file_data = new solar_mle_setup(trait_list, phenotype_filename,interp,  true); 
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
    Eigen::MatrixXd trait_matrix = file_data->return_output_matrix();
    Eigen::VectorXd trait_one =  trait_matrix.col(0);
    Eigen::VectorXd trait_two =  trait_matrix.col(1);
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
    Eigen::VectorXd final_parameters(8);
    double main_best_h;
    int main_iteration_count;
    double main_loglik = calculate_bivariate_model(trait_one,  trait_two,  S_ones, eigenvalues,  final_parameters,  main_best_h, main_iteration_count, 0,  debug); 
    if(main_loglik != main_loglik ){
        const char * error = "Convergence of this model could not be achieved";
        delete file_data;
        return error;
    }   
    if(display_pvalues){
        Eigen::VectorXd rhog_parameters(8);
        int iteration_count;
        double rhog_best_h;
        double rhog_loglik  = calculate_bivariate_model(trait_one,  trait_two,  S_ones, eigenvalues,  rhog_parameters,  rhog_best_h, iteration_count, 1,  debug);
        
        
        
        Eigen::VectorXd rhoe_parameters(8);
        double rhoe_best_h;
        double rhoe_loglik  = calculate_bivariate_model(trait_one,  trait_two,  S_ones, eigenvalues,  rhoe_parameters,  rhoe_best_h, iteration_count, 2,  debug);

       
        double rhog_pvalue = 1.0;
        if(rhog_loglik == rhog_loglik){
            rhog_pvalue = 2.0*chicdf(2.0*(main_loglik - rhog_loglik), 1); 
        }
        double rhoe_pvalue = 1.0;
        if(rhoe_loglik == rhoe_loglik){
            rhoe_pvalue = 2.0*chicdf(2.0*(main_loglik - rhoe_loglik), 1);
        }
        
        cout << "Final parameters\n";
        cout << "Pedigree: " << currentPed->filename() << endl;
        cout << "Phenotype: " << phenotype_filename << endl;
        cout << "Number of Subjects: " << eigenvalues.rows() << endl; 
        cout << "Total iterations: " << main_iteration_count << endl;
        cout << "Numerical Differentiation delta: " << main_best_h << endl;
        cout << trait_list[0] <<  "-mean: " << final_parameters(0) << " " << trait_list[1] << "-mean: "  << final_parameters(1) << endl;
        cout << trait_list[0] << "-h2r: " << calculate_constraint(final_parameters(2)) << " " << trait_list[1] << "-h2r: " << calculate_constraint(final_parameters(3)) << endl;
        cout << trait_list[0] << "-sd: "  << final_parameters(4) << " " << trait_list[1] << "-sd: " << final_parameters(5) << endl; 
        cout << "rhog: " << calculate_rho(final_parameters(6)) << " rhoe: " << calculate_rho(final_parameters(7))\
             << " rhop: " << calculate_rho(final_parameters(6))*sqrt(calculate_constraint(final_parameters(2)))*\
             sqrt(calculate_constraint(final_parameters(3))) + calculate_rho(final_parameters(7))*\
             sqrt(1.0-calculate_constraint(final_parameters(2)))*sqrt(1.0-calculate_constraint(final_parameters(3))) << endl;
        cout << "rhog p-value: " << rhog_pvalue <<  " rhog loglik: " << rhog_loglik << endl;
        cout << "rhoe p-value: " << rhoe_pvalue <<" rhoe loglik: " << rhoe_loglik << endl;
        cout << "loglik: " << main_loglik << endl;
       // Eigen::VectorXd errors = calculate_errors(trait_one, trait_two, S_ones,eigenvalues, final_parameters(0), final_parameters(1),\
                                        final_parameters(2), final_parameters(3), final_parameters(4), final_parameters(5), final_parameters(6), final_parameters(7), main_best_h);
        cout << endl;                                
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
        cout << "Final parameters\n";
        cout << "Pedigree: " << currentPed->filename() << endl;
        cout << "Phenotype: " << phenotype_filename << endl;
        cout << "Number of Subjects: " << eigenvalues.rows() << endl; 
        cout << "Total iterations: " << main_iteration_count << endl;
        cout << "Numerical Differentiation delta: " << main_best_h << endl;
        cout << trait_list[0] <<  "-mean: " << final_parameters(0) << " " << trait_list[1] << "-mean: "  << final_parameters(1) << endl;
        cout << trait_list[0] << "-h2r: " << calculate_constraint(final_parameters(2)) << " " << trait_list[1] << "-h2r: " << calculate_constraint(final_parameters(3)) << endl;
        cout << trait_list[0] << "-sd: "  << final_parameters(4) << " " << trait_list[1] << "-sd: " << final_parameters(5) << endl;
        cout << "rhog: " << calculate_rho(final_parameters(6)) << " rhoe: " << calculate_rho(final_parameters(7))\
             << " rhop: " << calculate_rho(final_parameters(6))*sqrt(calculate_constraint(final_parameters(2)))*\
             sqrt(calculate_constraint(final_parameters(3))) + calculate_rho(final_parameters(7))*\
             sqrt(1.0-calculate_constraint(final_parameters(2)))*sqrt(1.0-calculate_constraint(final_parameters(3))) << endl; 
        cout << "loglik: " << main_loglik << endl; 

        cout << endl;
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

    const char * error_msg = calculate_genetic_correlation(interp, trait_list,  phenotype_filename, get_pvalues, debug_mode, evd_data_filename);
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
            
            
