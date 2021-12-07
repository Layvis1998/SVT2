#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "Include/umfpack.h"
using namespace std;

double A = 20;    // mesh size
int n = 16;       // mesh of (n + 1) * (n + 1) nodes
int an = n + 1; 
int bn = n - 1;
int bbn = n - 2;
double h = A / n;
double sh = h * h;
int N = n * n;  
int aN = N + 1;
int a_N = an * an;

double dx = 0.01;
double dy = 1;
double porosity = 0.2;
double Pe = h / dx;

int* row_index = new int[aN];
int NNZ = 5 * bbn*bbn + 4 *(4 * bbn) + 4 * (3);                       
int* column_index = new int[NNZ];
double* values = new double[NNZ];
double* f_node = new double[a_N];
double* b = new double[N];
double *C_k = new double[N];
double *C_bk = new double[N];
double delta_t;
int cur = 0;

double* G_Top = new double[n];
double* G_Right = new double[n];
double* G_Bottom = new double[n];
double* G_Left = new double[n];

void No_borders(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx - h;
  cur++;

  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;

  values[cur] = 2 * (dx + dy) + sh / delta_t + h;
  column_index[cur] = i * n + j;
  cur++;
      
  values[cur] = - dy;
  column_index[cur] = i * n + (j + 1);
  cur++;
            
  values[cur] = - dx + h;
  column_index[cur] = (i + 1) * n + j ;
  cur++;
      
  row_index[k + 1] = row_index[k] + 5;
}

void Top_Dirichlet(int i, int j, int k)
{
  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
      
  values[cur] = 2 * (2 * dx + dy) + sh / delta_t + h;
  column_index[cur] = i * n + j;
  cur++;
      
  values[cur] = - dy;
  column_index[cur] = i * n + (j + 1);
  cur++;
            
  values[cur] = - dx;
  column_index[cur] = (i + 1) * n + j ;
  cur++;
      
  row_index[k + 1] = row_index[k] + 4; 
}

void Right_Dirichlet(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx - h;
  cur++;

  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;

  values[cur] = 2 * (dx + 2 * dy) + sh / delta_t + h;
  column_index[cur] = i * n + j;
  cur++;
            
  values[cur] = - dx;
  column_index[cur] = (i + 1) * n + j ;
  cur++;
      
  row_index[k + 1] = row_index[k] + 4;
}

void Left_Dirichlet(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx - h;
  cur++;

  values[cur] = 2 * (dx + 2 * dy) + sh / delta_t + h;
  column_index[cur] = i * n + j;
  cur++;
      
  values[cur] = - dy;
  column_index[cur] = i * n + (j + 1);
  cur++;
            
  values[cur] = - dx;
  column_index[cur] = (i + 1) * n + j;
  cur++;
      
  row_index[k + 1] = row_index[k] + 4;
}

void Bottom_Dirichlet(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx - h;
  cur++;

  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
      
  values[cur] = 2 * (2 * dx + dy) + sh / delta_t + h;
  column_index[cur] = i * n + j;
  cur++;
      
  values[cur] = - dy ;
  column_index[cur] = i * n + (j + 1);
  cur++;
    
  row_index[k + 1] = row_index[k] + 4;
}


void BR_Dirichlet(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx - h;
  cur++;

  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
      
  values[cur] = 4 * (dx + dy) + sh / delta_t + h;
  column_index[cur] = i * n + j;
  cur++;
      
  row_index[k + 1] = row_index[k] + 3;
}

void TL_Dirichlet(int i, int j, int k)
{
  values[cur] = 4 * (dx + dy) + sh / delta_t + h;
  column_index[cur] = i * n + j;
  cur++;      
      
  values[cur] = - dy;
  column_index[cur] = i * n + (j + 1);
  cur++;
            
  values[cur] = - dx;
  column_index[cur] = (i + 1) * n + j;
  cur++;
      
  row_index[k + 1] = row_index[k] + 3;
}

void TR_Dirichlet(int i, int j, int k)
{
  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
      
  values[cur] = 4 * (dx + dy) + sh / delta_t + h;
  column_index[cur] = i * n + j;
  cur++;
            
  values[cur] = - dx;
  column_index[cur] = (i + 1) * n + j;
  cur++;
      
  row_index[k + 1] = row_index[k] + 3;
}

void BL_Dirichlet(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx - h;
  cur++;

  values[cur] = 4 * (dx + dy) + sh / delta_t + h;
  column_index[cur] = i * n + j;
  cur++;
      
  values[cur] = - dy;
  column_index[cur] = i * n + (j + 1);
  cur++;
      
  row_index[k + 1] = row_index[k] + 3;  
}

double f(int i, int j)
{
  return 0.5;
}

void Fill_f_node()
{
  for (int a = 0; a < a_N; a++)
  {
    int i = a / an; 
    int j = a % an;
    f_node[i * an + j] = f(i, j);
  }
}

void Fill_b_from_f()
{
  for (int a = 0; a < N; a++)
  {
    int i = a / n; 
    int j = a % n;
    b[i * n + j] = f_node[i * an + j] + f_node[i * an + (j + 1)]
      + f_node[(i + 1) * an + j] + f_node[(i + 1) * an + (j + 1)];
    
    b[i * n + j] = b[i * n + j] / 4 * sh;
  }
}

void Fill_borders()
{
  for (int a = 0; a < n; a++)
    G_Left[a] = 0;
    
  for (int a = 0; a < n; a++)
    G_Right[a] = 0;  
  
  for (int a = 0; a < n; a++)
    G_Bottom[a] = 0;
    
  for (int a = 0; a < n; a++)
    G_Top[a] = 1;   
}

void Fill_b_from_borders()
{
  for (int a = 0; a < n; a++)
    b[a * n] += 2 * G_Left[a];
   
  for (int a = 0; a < n; a++)
    b[a * n + bn] += 2 * G_Right[a];  
  
  for (int a = 0; a < n; a++)
  {  
    b[bn * n + a] += 2 * G_Bottom[a];
    //b[bn * n + a] -= h / 2 * G_Bottom[a];
  }



  for (int a = 0; a < floor(n / 3 ); a++)
  {  
    b[a] += 2 * G_Top[a];       
  }

  for (int a = floor(n / 3 ); a < floor(2 * n / 3 ); a++)
  {  
    b[a] += 2 * G_Top[a];       
    b[a] += h * G_Top[a];  
  }
  
  for (int a = floor(2 * n / 3 ); a < n; a++)
  {  
    b[a] += 2 * G_Top[a];        
  }


}

void Print_C_k()
{
  cout << "Solution:\n ";
  for (int h = 0; h < n; h++)
  {
    uint32_t hght = h * n;
    cout << "[";
    for (int w = 0; w < bn; w++)
      cout << C_k[hght + w] << ", ";
    
    cout << C_k[hght + bn] << " ";
    cout<< "], \n";
  }

}

void Print_C_bk()
{
  cout << "Solution:\n ";
  for (int h = 0; h < n; h++)
  {
    uint32_t hght = h * n;
    cout << "[";
    for (int w = 0; w < bn; w++)
      cout << C_bk[hght + w] << ", ";
    
    cout << C_bk[hght + bn] << " ";
    cout<< "], \n";
  }

}

void Fill_b_from_C_bk()
{
  for (int a = 0; a < N; a++)
  {
    int i = a / n; 
    int j = a % n;
    b[i * n + j] += C_bk[i * n + j] / delta_t * sh;
  }
}


void Solve()
{
  //Factorizing  
  int start_fact = clock();

  void *Symbolic, *Numeric ;
  umfpack_di_symbolic (N, N, row_index, column_index, values, &Symbolic,
    nullptr, nullptr);
  umfpack_di_numeric (row_index, column_index, values, Symbolic, &Numeric,
    nullptr, nullptr);
  double time_fact = (double) (clock() - start_fact) / CLOCKS_PER_SEC;
  cout << "time of factorization: ms " << time_fact * 1e3 << endl;

  // Solving
  int start_solve = clock();
  umfpack_di_solve (UMFPACK_A, row_index, column_index, values, C_k, b,
    Numeric, nullptr, nullptr);
  double time_solve = (double )(clock() - start_solve) / CLOCKS_PER_SEC;

  umfpack_di_free_symbolic (&Symbolic);
  umfpack_di_free_numeric (&Numeric);  
  cout << "time of solving: ms " << time_solve * 1e3 << "\n\n";
}

int main(int argc, char** argv)
{
  cout << "Peclet number (it is equal to h/dx) = " << Pe << endl;

  int steps = atoi(argv[1]);
  delta_t = strtod(argv[2], NULL);
  for (int a = 0; a < N; a++)
    C_k[a] = 0;


  for (int i = 0; i < steps; i++)
  {
    row_index[0] = 0;
    cur = 0;  
    for (int j = 0; j < N; j++)
      C_bk[j] = C_k[j];

    Fill_f_node();
    Fill_b_from_f();
    Fill_borders();
    Fill_b_from_borders(); 
    Fill_b_from_C_bk();

    for (int k = 0; k < N; k++)
    {
      int i = k / n;
      int j = k % n;
    
      if ( (i < bn) && (j < bn) && (i > 0) && (j > 0) ) 
        No_borders(i, j, k);
    
      if ( (i == 0) && (j >= 1) && (j < bn) )
        Top_Dirichlet(i, j, k);
    
      if ( (j == bn) && (i >= 1) && (i < bn) )
        Right_Dirichlet(i, j, k);   
    
      if ( (j == 0) && (i >= 1) && (i < bn) ) 
        Left_Dirichlet(i, j, k);
    
      if ( (i == bn) && (j >= 1) && (j < bn) )
        Bottom_Dirichlet(i, j, k);

      if ( (i == bn) && (j == bn) ) //bottom right corner
        BR_Dirichlet(i, j, k);
        
      if ( (i == 0) && (j == 0) )   //top left corner      
        TL_Dirichlet(i, j, k);
     
      if ( (i == 0) && (j == bn) )  //top right corner
        TR_Dirichlet(i, j, k);
    
      if ( (i == bn) && (j == 0) )  //bottom left corner
        BL_Dirichlet(i, j, k);
    }

    Solve();
    Print_C_k();
  }  

}
