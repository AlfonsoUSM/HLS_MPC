// VARIABLE[ROWS][COLUMNS]
// N ROWS
// M COLUMNS

/* ADMM Algorithm
    // X Minimization
    v_x = z - c + u;
    x = R_inv * (-A.transpose()*v_x*rho - q);

    // Z Minimization
    z = (-A*x - u + c).max(0);

    // Update the scaled dual variable
    u = u + (A*x + z - c);
*/

/*
Variables:
    z
    u
    x
    -u
*/

/*
Parameters to save (fixed over iterations):
    -c
    R_inv
    -A^T*rho
    -q
    -A
    c
    A
*/

template<int N, int M, int P, typename T>
void mult(T (&A)[N][M], T (&B)[M][P], T (&res)[N][P])
{
	for(int i = 0; i < N; ++i)
	{
		for(int j = 0; j < P; ++j)
		{
			res[i][j] = A[i][0] * B[0][j];

			for(int k = 1; k < M; ++k)
			{
				res[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

template<int N, int M, typename T>
void add2(T (&A)[N][M], T (&B)[N][M], T (&res)[N][M])
{
	for(int i = 0; i < N; ++i)
	{
		for(int j = 0; j < M; ++j)
		{
			res[i][j] = A[i][j] + B[i][j];
		}
	}
}

template<int N, int M, typename T>
void add3(T (&A)[N][M], T (&B)[N][M], T (&C)[N][M], T (&res)[N][M])
{
	for(int i = 0; i < N; ++i)
	{
		for(int j = 0; j < M; ++j)
		{
			res[i][j] = A[i][j] + B[i][j] + C[i][j];
		}
	}
}

template<int N, int M, typename T>
void add4(T (&A)[N][M], T (&B)[N][M], T (&C)[N][M], T (&D)[N][M], T (&res)[N][M])
{
	for(int i = 0; i < N; ++i)
	{
		for(int j = 0; j < M; ++j)
		{
			res[i][j] = A[i][j] + B[i][j] + C[i][j] + D[i][j];
		}
	}
}

template<int N, int M, typename T>
void neg(T (&A)[N][M], T (&res)[N][M])
{
	for(int i = 0; i < N; ++i)
	{
		for(int j = 0; j < M; ++j)
		{
			res[i][j] = -A[i][j];
		}
	}
}

template<int N, int M, typename T>
void max(T (&A)[N][M], T sat_max)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (A[i][j] < sat_max)
            {
                A[i][j] = sat_max;
            }
        }
    }
}
//(T rho, T (&q_neg)[N][P], T (&rho_neg_A_T)[M][N], T (&c)[M][P], T (&R_inv)[N][N],
//T (&x)[N][P], T (&z)[M][P], T (&u)[M][P])
template<int N, int M, int P, typename T>
void x_min(T (&q_neg)[N][P], T (&rho_neg_A_T)[N][M], T (&c_neg)[M][P], T (&R_inv)[N][N],
            T (&x)[N][P], T (&z)[M][P], T (&u)[M][P])
{
    // v_x = z - c + u;
    T v_x[M][P];
    add3(z, c_neg, u, v_x);

    // x = R_inv * (-rho * A^T * v_x - q);
    T x1[N][P], x2[N][P], x3[N][P];
    mult(rho_neg_A_T, v_x, x1);
    add2(x1, q_neg, x2);
    mult(R_inv, x2, x);
}

template<int N, int M, int P, typename T>
void z_min(T (&A_neg)[M][N], T (&c)[M][P],
            T (&x)[N][P], T (&z)[M][P], T (&u)[M][P])
{
    T ub = 0;
    //  z = (-A*x - u + c).max(0);
    T loc_val1[M][P], u_neg[M][P];
    mult(A_neg, x, loc_val1);
    neg(u, u_neg);
    add3(loc_val1, u_neg, c, z);
    max(z, ub);
}

template<int N, int M, int P, typename T>
void u_update(T (&A)[M][N], T (&c_neg)[M][P],
            T (&x)[N][P], T (&z)[M][P], T (&u)[M][P])
{
    // u = u + (A*x + z - c);
    T loc_val1[M][P];
    mult(A,x,loc_val1);
    add4(u,loc_val1,z,c_neg,u);
}
