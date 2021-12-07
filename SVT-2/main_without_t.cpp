#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "Include/umfpack.h"
using namespace std;

int n = 25;        // mesh of (n + 1) * (n + 1) size
int an = n + 1; 
int bn = n - 1;
int bbn = n - 2;
double h = 1.0 / n;
double sh = h * h;
int N = n * n;  
int aN = N + 1;
int a_N = an * an;

int* row_index = new int[aN];
int NNZ = 5 * bbn*bbn + 4 *(4 * bbn) + 4 * (3);                       
int* column_index = new int[NNZ];
double* values = new double[NNZ];
double* b = new double[N];
int cur = 0;
  
int dx = 1;
int dy = 1;

void No_borders(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx;
  cur++;

  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;

  values[cur] = 2 * (dx + dy);
  column_index[cur] = i * n + j;
  cur++;
      
  values[cur] = - dy;
  column_index[cur] = i * n + (j + 1);
  cur++;
            
  values[cur] = - dx;
  column_index[cur] = (i + 1) * n + j;
  cur++;
      
  row_index[k + 1] = row_index[k] + 5;
}

void Top_Dirichlet(int i, int j, int k)
{
  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
      
  values[cur] = 2 * (2 * dx + dy);
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

void Right_Dirichlet(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx;
  cur++;

  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;

  values[cur] = 2 * (dx + 2 * dy);
  column_index[cur] = i * n + j;
  cur++;
            
  values[cur] = - dx;
  column_index[cur] = (i + 1) * n + j;
  cur++;
      
  row_index[k + 1] = row_index[k] + 4;
}

void Left_Dirichlet(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx;
  cur++;

  values[cur] = 2 * (2 * dx + dy);
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
  values[cur] = - dx;
  cur++;

  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
      
  values[cur] = 2 * (dx + 2 * dy);
  column_index[cur] = i * n + j;
  cur++;
      
  values[cur] = - dy;
  column_index[cur] = i * n + (j + 1);
  cur++;
    
  row_index[k + 1] = row_index[k] + 4;
}

void Top_Neumann(int i, int j, int k)
{
  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
      
  values[cur] = 2 * dx + dy;
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

void Right_Neumann(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx;
  cur++;

  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;

  values[cur] = dx + 2 * dy;
  column_index[cur] = i * n + j;
  cur++;
            
  values[cur] = - dx;
  column_index[cur] = (i + 1) * n + j;
  cur++;
      
  row_index[k + 1] = row_index[k] + 4;
}

void Left_Neumann(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx;
  cur++;

  values[cur] = dx + 2 * dy;
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

void Bottom_Neumann(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx;
  cur++;

  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
      
  values[cur] = dx + 2 * dy;
  column_index[cur] = i * n + j;
  cur++;
      
  values[cur] = - dy;
  column_index[cur] = i * n + (j + 1);
  cur++;
    
  row_index[k + 1] = row_index[k] + 4;
}

void BR_Dirichlet(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx;
  cur++;

  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
      
  values[cur] = 4 * (dx + dy);
  column_index[cur] = i * n + j;
  cur++;
      
  row_index[k + 1] = row_index[k] + 3;
}

void TL_Dirichlet(int i, int j, int k)
{
  values[cur] = 4 * (dx + dy);
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
      
  values[cur] = 4 * (dx + dy);
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
  values[cur] = - dx;
  cur++;

  values[cur] = 4 * (dx + dy);
  column_index[cur] = i * n + j;
  cur++;
      
  values[cur] = - dy;
  column_index[cur] = i * n + (j + 1);
  cur++;
      
  row_index[k + 1] = row_index[k] + 3;  
}

void TL_Neumann(int i, int j, int k)
{
  values[cur] = (dx + dy);
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

void TR_Neumann(int i, int j, int k)
{
  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
      
  values[cur] = (dx + dy);
  column_index[cur] = i * n + j;
  cur++;
            
  values[cur] = - dx;
  column_index[cur] = (i + 1) * n + j;
  cur++;
      
  row_index[k + 1] = row_index[k] + 3;
}

void BL_Neumann(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx;
  cur++;

  values[cur] = (dx + dy);
  column_index[cur] = i * n + j;
  cur++;
      
  values[cur] = - dy;
  column_index[cur] = i * n + (j + 1);
  cur++;
      
  row_index[k + 1] = row_index[k] + 3;  
}

void T_Dirichlet_L_Neumann(int i, int j, int k)
{
  values[cur] = 4 * dx + dy;
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

void T_Dirichlet_R_Neumann(int i, int j, int k)
{
  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
  
  values[cur] = 4 * dx + dy;
  column_index[cur] = i * n + j;
  cur++;
              
  values[cur] = - dx;
  column_index[cur] = (i + 1) * n + j;
  cur++;
      
  row_index[k + 1] = row_index[k] + 3;
}

void B_Dirichlet_R_Neumann(int i, int j, int k)
{
  values[cur] = - dx;
  column_index[cur] = (i - 1) * n + j;
  cur++;
  
  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
      
  values[cur] = 4 * dx + dy;
  column_index[cur] = i * n + j;
  cur++;
            
  row_index[k + 1] = row_index[k] + 3;
}

void B_Dirichlet_L_Neumann(int i, int j, int k)
{
  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
  
  values[cur] = 4 * dx + dy;
  column_index[cur] = i * n + j;
  cur++;
  
  column_index[cur] = (i + 1) * n + j;
  values[cur] = - dx;
  cur++;
      
  row_index[k + 1] = row_index[k] + 3;  
}

void B_Neumann_R_Dirichlet(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx;
  cur++;

  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
      
  values[cur] = dx + 4 * dy;
  column_index[cur] = i * n + j;
  cur++;
      
  row_index[k + 1] = row_index[k] + 3;  
}

void T_Neumann_L_Dirichlet(int i, int j, int k)
{
  values[cur] = dx + 4 * dy;
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

void T_Neumann_R_Dirichlet(int i, int j, int k)
{
  values[cur] = - dy;
  column_index[cur] = i * n + (j - 1);
  cur++;
      
  values[cur] = dx + 4 * dy;
  column_index[cur] = i * n + j;
  cur++;
            
  values[cur] = - dx;
  column_index[cur] = (i + 1) * n + j;
  cur++;
      
  row_index[k + 1] = row_index[k] + 3;
}

void B_Neumann_L_Dirichlet(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n + j;
  values[cur] = - dx;
  cur++;

  values[cur] = dx + 4 * dy;
  column_index[cur] = i * n + j;
  cur++;
      
  values[cur] = - dy;
  column_index[cur] = i * n + (j + 1);
  cur++;
      
  row_index[k + 1] = row_index[k] + 3;  
}

int main()
{
  double* G_Top = new double[n];
  double* G_Right = new double[n];
  double* G_Bottom = new double[n];
  double* G_Left = new double[n];
  double* f = new double[a_N];

  for (int a = 0; a < a_N; a++)
  {
    int i = a / an; 
    int j = a % an;
    f[i * an + j] = 0.5;
  }
  
  for (int a = 0; a < n; a++)
    G_Left[a] = 1;
    
  for (int a = 0; a < n; a++)
    G_Right[a] = 0;  
  
  for (int a = 0; a < n; a++)
    G_Bottom[a] = 0;
    
  for (int a = 0; a < n; a++)
    G_Top[a] = 0;   
  
  for (int a = 0; a < N; a++)
  {
    int i = a / n; 
    int j = a % n;
    b[i * n + j] = f[i * an + j] + f[i * an + (j + 1)] + f[(i + 1) * an + j]
      + f[(i + 1) * an + (j + 1)];
    b[i * n + j] = b[i * n + j] / 4 * sh;
  }
   
  for (int a = 0; a < n; a++)
    b[a * n] += 2 * G_Left[a];
   
  for (int a = 0; a < n; a++)
    b[a * n + bn] += 2 * G_Right[a];  
  
  for (int a = 0; a < n; a++)
    b[bn * n + a] += 2 * G_Bottom[a];
    
  for (int a = 0; a < n; a++)
    b[a] += 2 * G_Top[a];       
  
  delete [] G_Top;
  delete [] G_Right;
  delete [] G_Bottom;
  delete [] G_Left;
  delete [] f;  
  
  row_index[0] = 0;  
  for (int k = 0; k < N; k++)
  {
    int i = k / n;
    int j = k % n;
    
    if ( (i < bn) && (j < bn) && (i > 0) && (j > 0) ) 
      No_borders(i, j, k);
    
    if ( (i == 0) && (j >= 1) && (j < bn) )
      Top_Neumann(i, j, k);
    
    if ( (j == bn) && (i >= 1) && (i < bn) )
      Right_Dirichlet(i, j, k);   
    
    if ( (j == 0) && (i >= 1) && (i < bn) ) 
      Left_Dirichlet(i, j, k);
    
    if ( (i == bn) && (j >= 1) && (j < bn) )
      Bottom_Neumann(i, j, k);

    if ( (i == bn) && (j == bn) )  // bottom right corner
      B_Neumann_R_Dirichlet(i, j, k);
        
    if ( (i == 0) && (j == 0) )  //top left corner      
      T_Neumann_L_Dirichlet(i, j, k);
     
    if ( (i == 0) && (j == bn) ) //top right corner
      T_Neumann_R_Dirichlet(i, j, k);
    
    if ( (i == bn) && (j == 0) )  // bottom left corner
      B_Neumann_L_Dirichlet(i, j, k);
  }
  
  // Writing to file
  ofstream file("CSR.txt");
  file << N << " " << NNZ << endl;
  for (int k = 0; k < NNZ; ++k)
    file << values[k] << " ";
    
  file << endl;
  for (int k = 0; k < NNZ; ++k)
    file << column_index[k] << " ";
    
  file << endl;
  for (int i = 0; i < aN; ++i)
    file << row_index[i] << " ";

  file << endl;
  for (int i = 0; i < N; ++i)
    file << b[i] << " ";
    
  file << endl;
  file.close();

  //Factorizing  
  int start_fact = clock();
  double *u = new double[N];
  for (int i = 0; i < N; ++i)
    u[i] = 0;
    
  void *Symbolic, *Numeric ;
  umfpack_di_symbolic (N, N, row_index, column_index, values, &Symbolic, nullptr, nullptr);
  umfpack_di_numeric (row_index, column_index, values, Symbolic, &Numeric, nullptr, nullptr);
  double time_fact = (double) (clock() - start_fact) / CLOCKS_PER_SEC;
  cout << "time of factorization: ms " << time_fact * 1e3 << endl;

  // Solving
  int start_solve = clock();
  umfpack_di_solve (UMFPACK_A, row_index, column_index, values, u, b, Numeric, nullptr, nullptr);
  double time_solve = (double )(clock() - start_solve) / CLOCKS_PER_SEC;

  umfpack_di_free_symbolic (&Symbolic);
  umfpack_di_free_numeric (&Numeric);  
  cout << "time of solving: ms " << time_solve * 1e3 << "\n\n";

  delete[] b;
  delete[] values;
  delete[] column_index;
  delete[] row_index;

  //Checking
  cout << "Solution:\n ";
  for (int h = 0; h < n; h++)
  {
    uint32_t hght = h * n;
    cout << "[";
    for (int w = 0; w < n - 1; w++)
      cout << u[hght + w] << ", ";
    
    cout << u[hght + n - 1] << " ";
    cout<< "], \n";
  }
  delete[] u;
}
