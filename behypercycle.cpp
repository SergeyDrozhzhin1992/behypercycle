#include <fstream>
#include <iostream>
#include <stdio.h>
#include <cmath> 
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <glpk.h>
#include <algorithm>
#include <cstdlib>
using namespace std;

const double EPS = 1.0E-2;

/* Solve the linear programming problem: find the increments of the interaction matrix and the fixed point
  *
  *    INPUT:
  *        B_constr - matrix of restrictions on bij elements
  *        A - interaction matrix
  *        x - coordinates of a fixed point
  *        sizeA - the size of the interaction matrix
  *        f - mean fitness value
  *
  *    OUTPUT:
  *        B2 - matrix of (sizeA + 1) * sizeA
  *             first n lines - increments of the elements of the interaction matrix
  *             last line - increments of the fixed point
 */
gsl_matrix *solve_lin_prog(gsl_matrix *B_constr, gsl_matrix *A, gsl_vector *x, int sizeA, float f)
{
	/* Find the vector v */
    gsl_vector *v = gsl_vector_alloc(sizeA);
    for(int i = 0; i < sizeA; i++)
    {
        double s = 0;
        for(int j = 0; j < sizeA; j++)
            s = s + gsl_matrix_get(A, i, j) * gsl_vector_get(x, j);
        gsl_vector_set(v, i, s); 
	}
	    
	/* Find the matrix C */
	gsl_matrix *C = gsl_matrix_alloc(sizeA, sizeA);
	for(int i = 0; i < sizeA; i++)
		for(int j = 0; j < sizeA; j++)
			if (i == j) gsl_matrix_set(C, i, i, gsl_matrix_get(A, i, i) * gsl_vector_get(x, i) + gsl_vector_get(v, i));
			else gsl_matrix_set(C, i, j, gsl_matrix_get(A, i, j) * gsl_vector_get(x, i)); 
			
	
	/* Find the inverse matrix to the matrix C */
    gsl_matrix *invC = gsl_matrix_alloc(sizeA, sizeA);
    int s_c;  
    gsl_permutation *p_c = gsl_permutation_alloc(sizeA);
    gsl_linalg_LU_decomp(C, p_c, &s_c);
    gsl_linalg_LU_invert(C, p_c, invC);
    gsl_matrix_free(C);
    gsl_permutation_free(p_c);
	
    /* Find the inverse matrix to the matrix A */
    gsl_matrix *invA = gsl_matrix_alloc(sizeA, sizeA);
    gsl_matrix   *A2 = gsl_matrix_alloc(sizeA, sizeA);
    gsl_matrix_memcpy(A2, A);
    int s_a;
    
    gsl_permutation *p_a = gsl_permutation_alloc(sizeA);
    gsl_linalg_LU_decomp(A2, p_a, &s_a);
    gsl_linalg_LU_invert(A2, p_a, invA);
    gsl_matrix_free(A2);
    gsl_permutation_free(p_a);
        
	/* Find the constant */
	double const_c = 0;
	for(int k = 0; k < sizeA; k++)
	    for(int j = 0; j < sizeA; j++)
	        for(int i = 0; i < sizeA; i++)
	            const_c = const_c + gsl_matrix_get(invA, i, j) * gsl_matrix_get(invC, k, i);

    /* Find the coefficients before the elements bij */
    gsl_matrix *B = gsl_matrix_alloc(sizeA, sizeA);
    gsl_matrix_set_zero(B);

    for(int i = 0; i < sizeA; i++)
    {
        for(int j = 0; j < sizeA; j++)
        {

            double s1 = 0, s2 = 0;
            for(int l = 0; l < sizeA; l++)
            {
				s2 = s2 + gsl_matrix_get(invC, l, i) * gsl_vector_get(x, i) * gsl_vector_get(x, j);
                for(int m = 0; m < sizeA; m++)
                    for(int k = 0; k < sizeA; k++)
                        s1 = s1 + gsl_matrix_get(invA, k, i) * gsl_matrix_get(A, j, m) * gsl_vector_get(x, j) * gsl_matrix_get(invC, l, k) * gsl_vector_get(x, m);
			}
            gsl_matrix_set(B, i, j, (s1 + s2) / const_c);      

        }
    }  

    /* We set the linear programming problem */
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);

    /* Restrictions on the sum of bij and on the sum of aij * bij
     * In addition, if there are fixed point coordinates tending to 0 or 1,
     * impose restrictions on the increments of the corresponding coordinates
    */
    int count_chng = 0, count_chng2 = 2;
    for(int i = 0; i < sizeA; i++)
        if((gsl_vector_get(x, i) <= EPS) || (gsl_vector_get(x, i) >= (1 - EPS))) count_chng++;
    
    glp_add_rows(lp, 1 + count_chng);
    glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 0.0);
    
    if(count_chng > 0)
    { 
		for(int i = 0; i < sizeA; i++)
		{
			if(gsl_vector_get(x, i) <= EPS)
			{
				glp_set_row_bnds(lp, count_chng2, GLP_LO, 0, 0);
				count_chng2++;
			}
			
			if(gsl_vector_get(x, i) >= (1 - EPS))
			{
				glp_set_row_bnds(lp, count_chng2, GLP_UP, 0, 0);
				count_chng2++;
			}
			
			if((count_chng2 - 2) >= count_chng) break;
		}
    }

    int ia[(count_chng + 1) * sizeA * sizeA + 1], ja[(count_chng + 1) * sizeA * sizeA + 1];
    double ar[(count_chng + 1) * sizeA * sizeA + 1];
    int ind1 = 1;
    count_chng2 = 1;

    for(int k = 1; k <= (sizeA + 1); k++)
    {
		if((k == 1) || (gsl_vector_get(x, k - 2) <= EPS) || (gsl_vector_get(x, k - 2) >= (1 - EPS)))
		{
			for(int i = 0; i < sizeA; i++)
			{
				for(int j = 0; j < sizeA; j++)
				{
					ia[ind1] = count_chng2;
					ja[ind1] = i * sizeA + j + 1;
					
					if (count_chng2 == 1) 
						ar[ind1] = gsl_matrix_get(A, i, j);
					else
					{
						double s3 = 0, s4 = 0, s5 = 0;
						for(int l = 0; l < sizeA; l++)
						{
							for(int m = 0; m < sizeA; m++)
							{
								s3 = s3 + gsl_matrix_get(invC, k - 2, m) * gsl_matrix_get(invA, m, l);
								s4 = s4 + gsl_matrix_get(invC, k - 2, m) * gsl_matrix_get(invA, m, l) * gsl_matrix_get(A, l, i) * gsl_vector_get(x, i) * gsl_vector_get(x, j);
								s5 = s5 + gsl_matrix_get(invC, k - 2, m) * gsl_matrix_get(invA, m, i) * gsl_matrix_get(A, j, l) * gsl_vector_get(x, l) * gsl_vector_get(x, j);
							}
						}
						ar[ind1] = s3 * gsl_matrix_get(B, i, j) - s4 - s5;
					}
					
					ind1 = ind1 + 1;        
				}
			}
			count_chng2++;
	    }
	}

    /* Restrictions on each bij left and right */
    ind1 = 1;
    glp_add_cols(lp, sizeA * sizeA);
    for(int i = 0; i < sizeA; i++)
    {
        for(int j = 0; j < sizeA; j++)
        {
            glp_set_col_bnds(lp, ind1, GLP_DB, -gsl_matrix_get(B_constr, i, j), gsl_matrix_get(B_constr, i, j));
            ind1 = ind1 + 1;        
        }
    }

    /* The coefficients in the equation that we want to maximize */
    ind1 = 1;
    for(int i = 0; i < sizeA; i++)
    {
        for(int j = 0; j < sizeA; j++)
        {
            glp_set_obj_coef(lp, ind1, gsl_matrix_get(B, i, j));
            ind1 = ind1 + 1;        
        }
    }

    /* Solve the linear programming problem */
    glp_load_matrix(lp, (count_chng + 1) * sizeA * sizeA, ia, ja, ar);
    glp_simplex(lp, NULL);

    ind1 = 1;
    for(int i = 0; i < sizeA; i++)
    {
        for(int j = 0; j < sizeA; j++)
        {
            gsl_matrix_set(B, i, j, glp_get_col_prim(lp, ind1));
            ind1 = ind1 + 1;        
        }
    }
    
    /*Find the increment of the fixed point*/
    gsl_vector *u = gsl_vector_alloc(sizeA);
    double z = glp_get_obj_val(lp);
    
    for(int k = 0; k < sizeA; k++)
    {
		double s6 = 0, s7 = 0;
        for(int l = 0; l < sizeA; l++)
            for(int j = 0; j < sizeA; j++)
            {
				s7 = s7 + gsl_matrix_get(invC, k, j) * gsl_matrix_get(invA, j, l);
                for(int i = 0; i < sizeA; i++)
                   for(int m = 0; m < sizeA; m++)
                       s6 = s6 + gsl_matrix_get(invC, k, m) * gsl_matrix_get(invA, m, l) * (gsl_matrix_get(A, l, i) * gsl_vector_get(x, i) * gsl_matrix_get(B, i, j) + gsl_matrix_get(B, l, i) * gsl_vector_get(x, i) * gsl_matrix_get(A, i, j)) * gsl_vector_get(x, j);
	        }
	    gsl_vector_set(u, k, z * s7 - s6);
	}
	gsl_matrix_free(invA);
	gsl_matrix_free(invC);
	gsl_vector_free(v);

    glp_delete_prob(lp);
    
    gsl_matrix *B2 = gsl_matrix_alloc(sizeA + 1, sizeA);
    for(int i = 0; i < sizeA; i++)
    {
		for(int j = 0; j < sizeA; j++)
		    gsl_matrix_set(B2, i, j, gsl_matrix_get(B, i, j));
		gsl_matrix_set(B2, sizeA, i, gsl_vector_get(u, i));
	}  

    gsl_matrix_free(B);
    gsl_vector_free(u);
    return B2;
}

/* Solve ODE */
int func (double t, const double y[], double f[], void *params)
{
    (void)(t); 
    gsl_matrix *A = (gsl_matrix*)params;
        
    for(int k = 0; k < A->size1; k++)
    {
		f[k] = 0;
		double s1 = 0, s2 = 0;
        for(int m = 0; m < A->size1; m++)
        {
			for(int i = 0; i < A->size1; i++)
			{
				for(int j = 0; j < A->size1; j++)
					s2 = s2 + gsl_matrix_get(A, m, j) * gsl_matrix_get(A, j, i) * y[j] * y[i] * y[m];
				s1 = s1 + gsl_matrix_get(A, k, i) * gsl_matrix_get(A, i, m) * y[i] * y[m];
			}
		}
		f[k] = y[k] * (s1 - s2);
	}
    
    return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    (void)(t); 
    gsl_matrix *A = (gsl_matrix*)params;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, A->size1, A->size1);
    gsl_matrix * m = &dfdy_mat.matrix;
    
    for(int k = 0; k < A->size1; k++)
    {
		for(int l = 0; l < A->size1; l++)
		{
			
			double s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0;
			for(int i = 0; i < A->size1; i++)
			{
				for(int j = 0; j < A->size1; j++)
				{
					s3 = s3 + gsl_matrix_get(A, l, j) * gsl_matrix_get(A, j, i) * y[j] * y[i];
					s4 = s4 + gsl_matrix_get(A, i, j) * gsl_matrix_get(A, j, l) * y[j] * y[i];
					s5 = s5 + gsl_matrix_get(A, i, l) * gsl_matrix_get(A, l, j) * y[j] * y[i];
				}
				s1 = s1 + gsl_matrix_get(A, k, i) * gsl_matrix_get(A, i, l) * y[i];
				s2 = s2 + gsl_matrix_get(A, k, l) * gsl_matrix_get(A, l, i) * y[i];
			}
			
			double s6 = 0, s7 = 0;
			for(int g = 0; g < A->size1; g++)
			{
				for(int i = 0; i < A->size1; i++)
				{
					for(int j = 0; j < A->size1; j++)
						s6 = s6 + gsl_matrix_get(A, g, j) * gsl_matrix_get(A, j, i) * y[j] * y[i] * y[g];
					s7 = s7 + gsl_matrix_get(A, k, i) * gsl_matrix_get(A, i, g) * y[i] * y[g];
				}
			}
			
			if(k == l) gsl_matrix_set(m, k, l, s7 - s6 + y[k] * (s1 + s2 - s3 - s4 - s5));
			else gsl_matrix_set(m, k, l, y[k] * (s1 + s2 - s3 - s4 - s5));
		}
	}
    
    for(int i = 0; i < A->size1; i++)
        dfdt[i] = 0.0;
 
    return GSL_SUCCESS;
}

/* Calculate average integral fitness (quadrature formulas) */
double get_avg_integral_fitness(gsl_matrix *U_continuos, gsl_matrix *A, int sizeA, int count_solve_step2, double count_step)
{
	double s, f = 0;
	for(int i = 0; i <= count_step; i++)
	{
		s = 0;
		for(int j = 0; j < sizeA; j++)
		    for(int k = 0; k < sizeA; k++)
		        for(int m = 0; m < sizeA; m++)
		            s = s + gsl_matrix_get(A, j, m) * gsl_matrix_get(A, m, k) * gsl_matrix_get(U_continuos, count_solve_step2, j * (count_step + 1) + i) * gsl_matrix_get(U_continuos, count_solve_step2, k * (count_step + 1) + i) * gsl_matrix_get(U_continuos, count_solve_step2, m * (count_step + 1) + i);
		
		if ((i == 0) || (i == count_step)) s = s / 2;
		f = f + s;
    }
	   
	f = f / count_step;
	
	return f;
}


/* Write data to a file just to see it with your eyes:
  * write to the file the evolution of the interaction matrix, 
  * the matrix norm at different iterations 
  * and the fitness system
 */
void write_in_file(int sizeA, int count_iter, gsl_matrix *A_time, gsl_vector *matrix_norm_vec, gsl_vector *fitness_vec)
{
	
	/* Write the type of matrix A at each step */
	ofstream evolution_A("evolution_matrix_A.txt");
	for(int i = 0; i <= count_iter; i++)
	{
		for(int j = 0; j < sizeA; j++)
		{
			for(int k = 0; k < sizeA; k++)
			{
				evolution_A << gsl_matrix_get(A_time, j * sizeA + k, i) << " ";
			}
			evolution_A << endl;
		}
		evolution_A << endl;
	}
	evolution_A.close();
	
	/* Write the norm of the matrix A at each step */
	ofstream norm_A("norma_matrix_A.txt");
	for(int i = 0; i <= count_iter; i++)
		norm_A << gsl_vector_get(matrix_norm_vec, i) << endl;
	norm_A.close();
	
	/* Write fitness at every step */
	ofstream fitness("fitness.txt");
	for(int i = 0; i <= count_iter; i++)
		fitness << gsl_vector_get(fitness_vec, i) << endl;
	fitness.close();
}


/* Write the data to a file in a binary form, then we consider it as Matlab */
void write_in_file_for_Matlab(int sizeA, int count_iter, double count_step, int solve_step, double count_solve_step, 
                              gsl_matrix *U, gsl_vector *fitness_vec, gsl_vector *fitness_vec_avg, gsl_matrix *U_continuos, 
                              gsl_vector *time_vec, gsl_matrix *A_time)
{
	double num;
	
	/* Write the evolution of the matrix */
	ofstream evolution_A("evolution_matrix_matlab.txt", ios::binary | ios::out);
	
	for(int j = 0; j < sizeA; j++)
		for(int k = 0; k < sizeA; k++)
		    for(int i = 0; i <= count_iter; i++)
		    {
				num = gsl_matrix_get(A_time, j * sizeA + k, i); 
				evolution_A.write((char*)&num, sizeof num);
			}
	evolution_A.close();
	
	/* Write frequencies */
	ofstream freq("freqType_matlab.txt", ios::binary | ios::out);
	for(int i = 0; i < sizeA; i++)
	{
		for(int j = 0; j <= count_iter; j++)
		{
			num = gsl_matrix_get(U, i, j); 
			freq.write((char*)&num, sizeof num);
		}
	}
	freq.close();
	
	/* Write the fitness vector */
	ofstream fitn("fitness_matlab.txt", ios::binary | ios::out);
	for(int i = 0; i <= count_iter; i++)
	{
		num = gsl_vector_get(fitness_vec, i);
		fitn.write((char*)&num, sizeof num);
	}				
	fitn.close();
	
	/* Write the vector of average integral fitness */
	ofstream fitn_avg("fitness_avg_matlab.txt", ios::binary | ios::out);
	for(int i = 0; i < count_solve_step; i++)
	{
		num = gsl_vector_get(fitness_vec_avg, i);
		fitn_avg.write((char*)&num, sizeof num);
	}				
	fitn_avg.close();
	
    /* Write the decisions of the ODE */
	ofstream freq_cont("freqType_continuos_matlab.txt", ios::binary | ios::out);
	for(int i = 0; i < count_solve_step; i++)
	{
		for(int j = 0; j < sizeA * (count_step + 1); j++)
		{
				num = gsl_matrix_get(U_continuos, i, j);
				freq_cont.write((char*)&num, sizeof num);
		}
	}
	freq_cont.close();
	
	/* Write the time vector */
	ofstream time("time_matlab.txt", ios::binary | ios::out);
	for(int i = 0; i <= count_step; i++)
	{
		num = gsl_vector_get(time_vec, i);
		time.write((char*)&num, sizeof num);
	} 
	time.close();
	
	/* Write custom data for MatLab */
	ofstream set("settings_matlab.txt");
	set.write((char*)&sizeA, sizeof sizeA);
	set.write((char*)&count_iter, sizeof count_iter);
	set.write((char*)&count_step, sizeof count_step);
	set.write((char*)&solve_step, sizeof solve_step);
	set.write((char*)&count_solve_step, sizeof count_solve_step);
	set.close();
}


int main(int *argc, char **argv)
{
   /* Enter the data from the keyboard:
    *     interaction matrix size
    *     the number of iterations of evolution
    *     the final point in time to solve the ODE
    *     time step (grid for solving ODE)
    *     ODU decision step (at which iterations of evolution we will solve an ODE) 
   */
    int sizeA, count_iter, solve_step;
    double t1, h;
    /*interaction matrix size*/
    cout << "Enter the size of the matrix A "; cin >> sizeA; cout << endl;
    /*the number of iterations of evolution*/
    cout << "Enter the count iteration "; cin >> count_iter; cout << endl;
    /*the final point in time to solve the ODE*/
    cout << "Enter T1 "; cin >> t1; cout << endl;
    /*time step (grid for solving ODE)*/
    cout << "Enter the time step "; cin >> h; cout << endl;
    /*ODU decision step (at which iterations of evolution we will solve an ODE) */
    cout << "Enter the solve dif.eq. step "; cin >> solve_step; cout << endl;    

    /* We read from the files the interaction matrix, the restrictions on the bij elements, as well as the initial conditions */
    gsl_matrix *A = gsl_matrix_alloc(sizeA, sizeA);
    gsl_matrix *B_constr = gsl_matrix_alloc(sizeA, sizeA);
    gsl_vector *u0 = gsl_vector_alloc(sizeA);
    
    double buff;
    ifstream fin_A("Matrix_A.txt");// opened the file for reading the interaction matrix
    ifstream fin_B_constr("Matrix_B_constr.txt");// opened the file to read the restrictions on the elements bij
    ifstream fin_u0("u0.txt");// opened the file to read the initial data
    for(int i = 0; i < sizeA; i++)
    {
		fin_u0 >> buff; gsl_vector_set(u0, i, buff);
        for(int j = 0; j < sizeA; j++)
        {
            fin_A         >> buff;   gsl_matrix_set(A,               i, j,    buff);
            fin_B_constr  >> buff;   gsl_matrix_set(B_constr,        i, j,    buff);          
        }
    }        
    fin_A.close();   fin_B_constr.close(); fin_u0.close();// close the file
    

    /*Output*/
    //equilibrium matrix at each iteration
    gsl_matrix *U = gsl_matrix_alloc(sizeA, count_iter + 1);
    //average fitness vector at each iteration
    gsl_vector *fitness_vec = gsl_vector_alloc(count_iter + 1);
    //Matrix of the evolutionary change of the elements of the interaction matrix. 
    //The first column is the elements of the interaction matrix at the initial time (a11, a12, a1n, ..., an1, ..., ann), etc.
    gsl_matrix *A_time = gsl_matrix_alloc(sizeA * sizeA, count_iter + 1); 
    //matrix norm at each iteration
    gsl_vector *matrix_norm_vec = gsl_vector_alloc(count_iter + 1);
    
    
    /* Data for solving ODE */
    double y[sizeA], count_step, t0, count_solve_step;
    modf(t1 / h, &count_step);
    modf((count_iter + 1) / solve_step, &count_solve_step); 
    if((count_solve_step * solve_step == (count_iter + 1)) || (count_iter == 0)) count_solve_step++;
    else count_solve_step = count_solve_step + 2;
    int solve_step2 = 1, count_solve_step2 = 0;
    gsl_matrix *U_continuos = gsl_matrix_alloc(count_solve_step, sizeA * (count_step + 1));
    gsl_odeiv2_driver * d;
    gsl_odeiv2_system sys;
    
    /* Data for finding average integral fitness */
    gsl_vector *fitness_vec_avg = gsl_vector_alloc(count_solve_step);
    
    gsl_vector *time_vec = gsl_vector_alloc(count_step + 1); 
    for(int i = 0; i <= count_step; i++)
        gsl_vector_set(time_vec, i, h * i);
    
    gsl_vector *x = gsl_vector_alloc(sizeA);    
    gsl_matrix *B; 
    
    /* At each next step, we again find the frequencies and fitness, solve the LPP and get a new matrix A */
    for(int i = 0; i <= count_iter; i++)
    {
		cout << "I = " << i << endl << endl; /* Print the iteration number */
        
        /* The value of the fixed point at 0 iteration */
        if(i == 0) gsl_vector_set_all(x, 1.0 / sizeA);
        
        if(gsl_vector_min(x) >= 0)
        {
			gsl_matrix_set_col(U, i, x);  

            /* Calculate Fitness */
			gsl_vector_set(fitness_vec, i, 0);
			for(int k = 0; k < sizeA; k++)
				for(int j = 0; j < sizeA; j++)
					for(int m = 0; m < sizeA; m++)
					    gsl_vector_set(fitness_vec, i, gsl_vector_get(fitness_vec, i) + gsl_matrix_get(A, k, m) * gsl_matrix_get(A, m, j) * gsl_vector_get(x, k) * gsl_vector_get(x, j) * gsl_vector_get(x, m));					    				
					
			/* Calculate the matrix norm */
			gsl_vector_set(matrix_norm_vec, i, 0);
			for(int k = 0; k < sizeA; k++)
				for(int j = 0; j < sizeA; j++)
				    gsl_vector_set(matrix_norm_vec, i, gsl_vector_get(matrix_norm_vec, i) + gsl_matrix_get(A, k, j) * gsl_matrix_get(A, k, j));
	 
			/* Save the values of the elements of the interaction matrix */
			for(int k = 0; k < sizeA; k++)
				for(int j = 0; j < sizeA; j++)
					gsl_matrix_set(A_time, k * sizeA + j, i, gsl_matrix_get(A, k, j));
						
			/* Solve the ODE and write the result in the matrix */
			if (i == 0 || i == count_iter || solve_step2 == solve_step)
			{
				sys = {func, jac, sizeA, A};
				d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
			
				for(int k = 0; k < sizeA; k++) 
				{
					y[k] = gsl_vector_get(u0, k);
					gsl_matrix_set(U_continuos, count_solve_step2, k * (count_step + 1), y[k]);
				}
			
				t0 = 0.0;
				for(int k = 0; k < count_step; k++)
				{
					double ti = h * (k + 1);
					int status = gsl_odeiv2_driver_apply (d, &t0, ti, y);
					if (status != GSL_SUCCESS)
					{
						printf ("error, return value=%d\n", status);
						break;
					}
					for(int j = 0; j < sizeA; j++)
					    gsl_matrix_set(U_continuos, count_solve_step2, j * (count_step + 1) + k + 1, y[j]); 
				}
				
				/* Calculate average integral fitness */
				gsl_vector_set(fitness_vec_avg, count_solve_step2, get_avg_integral_fitness(U_continuos, A, sizeA, count_solve_step2, count_step));
				
				gsl_odeiv2_driver_free(d);  
				count_solve_step2++;
			}      
			solve_step2++;
			if (solve_step2 > solve_step) solve_step2 = 1;
	    }
	    
	    else
	    {
		    count_iter = i - 1;
		    count_solve_step = count_solve_step2;
		    cout << "COMPONENTS LESS 0" << endl;
		    for(int ii = 0; ii < sizeA; ii++)
		        cout << gsl_vector_get(x, ii) << "  ";
		    break;    	
		}
		    
		 
        /* We solve the linear programming problem: 
         *     we find the increments of the elements of the interaction matrix and the fixed point 
        */
        B = solve_lin_prog(B_constr, A, x, sizeA, gsl_vector_get(fitness_vec, i));
        
        gsl_matrix *B2 = gsl_matrix_alloc(sizeA, sizeA);
        for(int k = 0; k < sizeA; k++)
            for(int j = 0; j < sizeA; j++)
                gsl_matrix_set(B2, k, j, gsl_matrix_get(B, k, j));
                
        gsl_vector *x2 = gsl_vector_alloc(sizeA);
        for(int k = 0; k < sizeA; k++)
            gsl_vector_set(x2, k, gsl_matrix_get(B, sizeA, k));    
              
        /* Rewrite the matrix A */
        gsl_matrix_add(A, B2);
        gsl_vector_add(x, x2);
        
        gsl_vector_free(x2);
        gsl_matrix_free(B);
        gsl_matrix_free(B2);
    }   
    
    write_in_file(sizeA, count_iter, A_time, matrix_norm_vec, fitness_vec);
	write_in_file_for_Matlab(sizeA, count_iter, count_step, solve_step, count_solve_step, U, fitness_vec, fitness_vec_avg, U_continuos, time_vec, A_time); 

    gsl_matrix_free(A);
    gsl_vector_free(matrix_norm_vec);
    gsl_matrix_free(B_constr); 
    gsl_matrix_free(U);
    gsl_matrix_free(A_time);
    gsl_matrix_free(U_continuos);
    gsl_vector_free(fitness_vec);
    gsl_vector_free(fitness_vec_avg);
    gsl_vector_free(time_vec);
    gsl_vector_free(u0);
    gsl_vector_free(x);
            
    return 0;
}
