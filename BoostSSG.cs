using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace pssaclass
{
    class BoostSSG
    {
               int N;
        int L;
        int K;

        double[] A;
        double[] AT;
        double[] C;

        double[] S;
        double[] U;
        double[] VT;

        int CN;

        public BoostSSG(double[] Signal, int Len, int Window, double[] Spectrum, int ComponentsNum)
        {
            N = Len;
            L = Window;
            K = N-L+1;

            CN = ComponentsNum;

            A = new double[L * K];
            AT = new double[K * L];
            C = new double[L * L];

            S = new double[CN];
            U = new double[L * CN];
            VT = new double[L * CN];

            //Track matrix
            for (int j = 0; j < K; j++)
            {
                for (int i = 0; i < L; i++)
                {
                    A[i * K + j] = Signal[(j + i)];
                }
            }
/*
            Console.Write("A=\n");
            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    Console.Write("{0:f6}\t", A[i * K + j]);
                }
                Console.Write("\n");
            }
            Console.Write("\n");
*/

            //Step 1. Transpose A -> AT
            this.TransposeMatrixThisResultOfTransposeMatrix(A, L, K, AT);
            //Step 2. X=A*AT
            this.MatrixCThisResultOfMultiplyingMatrixAOnMatrixB(C, A, L, K, AT, K, L);
/*
            Console.Write("C=\n");
            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < L; j++)
                {
                    Console.Write("{0:f6}\t", C[i * L + j]);
                }
                Console.Write("\n");
            }
            Console.Write("\n");
*/


            this.Eigensrch(C, L, S, U);

            //this.MatrixCThisResultOfMultiplyingMatrixAOnMatrixB(VT, X, K, L, U, L, CN);

/*
            for (int j = 0; j < CN; j++)
            {
                for (int i = 0; i < K; i++)
                {
                    VT[i * CN + j] = VT[i * CN + j] / Math.Sqrt(Math.Abs(S[j]));
                }
            }

            Array.Clear(Spectrum, 0, Spectrum.Length);

            //for (int c = 0; c < CN; c++)
            Parallel.For(0, CN, c =>
            {
                double[] R = new double[K * K];

                for (int i = 0; i < K; i++)
                {
                    for (int j = 0; j < K; j++)
                    {
                        R[i * K + j] = U[i * CN + c] * Math.Sqrt(Math.Abs(S[c])) * VT[j * CN + c];
                    }
                }

                for (int j = 0; j < L; j++)
                {
                    for (int i = 0; i < (j + 1); i++)
                    {
                        Spectrum[N * c + j] = Spectrum[N * c + j] + R[i * K + (j - i)];
                    }
                    Spectrum[N * c + j] = Spectrum[N * c + j] / (j + 1);
                }

                for (int j = L; j < K; j++)
                {
                    for (int i = 0; i < L; i++)
                    {
                        Spectrum[N * c + j] = Spectrum[N * c + j] + R[i * K + (j - i)];
                    }
                    Spectrum[N * c + j] = Spectrum[N * c + j] / L;
                }

                for (int j = K; j < N; j++)
                {
                    for (int i = 0; i < L - 0 - (j - K); i++)
                    {
                        Spectrum[N * c + j] = Spectrum[N * c + j] + R[(j - K + i) * K + (K - 1 - i)];
                    }
                    Spectrum[N * c + j] = Spectrum[N * c + j] / (N - j);
                }

            });

 */
        }

        private void MatrixTAndMatrixQIsResultOfLanczosTridiagonalizationMatrixA(double[] A, int RS, int CS, double[] T, double[] Q)
        {
            double[] z = new double[RS];
            int m = 0;

            Array.Clear(z, 0, z.Length);
            Array.Clear(T, 0, T.Length);
            Array.Clear(Q, 0, Q.Length);

            Q[0] = 1.0;

            for (int i = 0; i < CN; i++)
            {
                for (int n = 0; n < RS; n++)
                {
                    z[n] = (double)0.0;
                    for (int j = 0; j < CS; j++)
                    {
                        z[n] = z[n] + A[n * CS + j] * Q[j * CN + i];
                    }
                }

                T[CS * i + i] = (double)0.0;

                for (int n = 0; n < RS; n++)
                {
                    T[CN * i + i] = T[CN * i + i] + Q[n * CN + i] * z[n];
                }

                if (i != CN - 1)
                {
                    for (int n = 0; n < RS; n++)
                    {
                        Q[n * CN + i + 1] = (i == 0) ? z[n] - T[CN * i + i] * Q[n * CN + i] : z[n] - T[CN * i + i] * Q[n * CN + i] - T[(i - 1) * CN + i] * Q[n * CN + i - 1];
                    }

                    double R;
                    for (int n = 0; n < m; n++)
                    {
                        R = 0;
                        for (int k = 0; k < RS; k++)
                        {
                            R = R + Q[k * CN + n] * Q[k * CN + i + 1];
                        }
                        for (int k = 0; k < RS; k++)
                        {
                            Q[k * CN + i + 1] = Q[k * CN + i + 1] - Q[k * CN + n] * R;
                        }
                    }

                    T[i * CN + i + 1] = 0.0;
                    for (int n = 0; n < RS; n++)
                        T[i * CN + i + 1] = T[i * CN + i + 1] + Q[n * CN + i + 1] * Q[n * CN + i + 1];
                    T[i * CN + i + 1] = Math.Sqrt(T[i * CN + i + 1]);

                    for (int n = 0; n < RS; n++)
                        Q[n * CN + i + 1] = Q[n * CN + i + 1] / T[i * CN + i + 1];

                    T[(i + 1) * CN + i] = T[i * CN + i + 1];
                }
                m++;
            }
            z = null;
        }

        private void TransposeMatrixThisResultOfTransposeMatrix(double[] Matrix, int L, int K, double[] TransposeMatrix)
        {
            for (int i = 0; i < L; i++)
            //Parallel.For(0, L, i =>
            {
                for (int j = 0; j < K; j++)
                {
                    TransposeMatrix[j * L + i] = Matrix[i * K + j];
                }
            }
        }

        private void MatrixCThisResultOfMultiplyingMatrixAOnMatrixB(double[] C, double[] A, int AN, int AM, double[] B, int BN, int BM)
        {
            Array.Clear(C, 0, C.Length);

            //for (int i = 0; i < AN; i++)
            Parallel.For(0, AN, i =>
            {
                for (int j = 0; j < BM; j++)
                {
                    for (int k = 0; k < BN; k++)
                    {
                        C[i * BM + j] = C[i * BM + j] + A[i * AM + k] * B[k * BM + j];
                    }
                }
            });
        }

        private void GetApBp(double[] T, double[] a, double[] b, int n)
        {
            for (int i = 0; i < n; i++)
            {
                b[i]=T[i * n  + i];
            }
            a[0] = 0;
            for (int i = 0; i < n-1; i++)
            {
                a[i+1] = T[i * n + i+1];
            }
        }

        private void GetC_Сoef(double[] c_coef, double[] a, double[] b, int n)
        {
            for (int i = 0; i < n; i++)
            {
                c_coef[i] = (i == 0) ? 0 : -(a[i] / (b[i - 1] + c_coef[i - 1] * a[i - 1]));
            }
        }

        private void GetQp(double[] q, double[] a, double[] b, double[] c_coef, int n)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == j) q[i * n + j] = 1.0;
                    else
                        if ((j < i) && (i <= n - 2))
                        {
                            double mull = 1.0;
                            for (int k = j + 2; k < i + 2; k++)
                            {
                                mull =  c_coef[k]*mull;
                            }
                            q[i * n + j] = (a[j + 1] / a[i + 1]) * mull;   
                        }
                        else
                            if ((i == n - 1) && (0 <= j) && (j <= n - 3))
                            {
                                double mull = 1.0;
                                for (int k = j + 2 ; k < n ; k++)
                                {
                                    mull = mull * c_coef[k];
                                }
                                q[i * n + j] = (-a[j + 1] / (b[n - 1] + c_coef[n - 1] * a[n - 1])) * mull;
                            }

                            else
                                if ((i == n - 1) && (j == n - 2)) q[i * n + j] = (-a[n - 1] / (b[n - 1] + c_coef[n - 1] * a[n - 1]));
                                else
                                    q[i * n + j] = 0.0;                     
                }
            }
        }

        private void GetCp(double[] c, double[] a, double[] b, double[] c_coef, int n)
        {
            for (int i = n-1; i > -1; i--)
            {              
                for (int j = n-1; j > -1; j--)
                {
                    if (i == j) c[i * n + j] = 1.0;
                    else                  
                        if (i < j)
                        {
                            double mull = 1.0;
                            for (int k = j; k > i; k--)
                            {
                                mull = c_coef[k]*mull;                             
                            }
                            c[i * n + j] = mull;
                        }
                        else
                            c[i * n + j] = 0.0;
                }
            }
        }

        private void GetDp(double[] d, double[] a, double[] b, double[] c_coef, int n)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if ((i == 0) && (j == 0)) d[i * n + j] = 1 / b[0];
                    else
                        if ((1 <= i) && (j <= n - 2) && (i == j))
                        {
                            d[i * n + j] = -c_coef[i + 1] / a[i + 1];
                        }
                        else
                            if ((i == n - 1) && (j == n - 1))
                                d[i * n + j] = 1.0 / (b[n - 1] + c_coef[n - 1] * a[n - 1]);
                            else
                                d[i * n + j] = 0.0;
                }
            }
        }

        private void MulMatrixOnVector(double[] A, int RS, int CS, double[] u1, double[] vec)
        {
            Array.Clear(vec, 0, vec.Length);

            for (int i = 0; i < RS; i++)
            {
                for (int j = 0; j < CS; j++)
                {
                    vec[i] = vec[i] + A[i * CS + j] * u1[j];
                }
            }
        }

        private double NormVector(double[] u1, int L)
        {
            double norm = 0;
            for (int i = 0; i < L; i++) norm = norm + u1[i] * u1[i];
            norm = Math.Sqrt(norm);
            for (int i = 0; i < L; i++) u1[i] = (double)(u1[i] / norm);

            return norm;
        }

        private double ComputeEigenU1Lambda1(double[] A, int RS, int CS, double[] u1, double epsilon)
        {
            double norm1;
            double norm2;
            double[] vec = new double[RS];
            norm2 = NormVector(u1, RS);
            do
            {
                norm1 = norm2;
                MulMatrixOnVector(A, RS, CS, u1, vec);
                Array.Copy(vec, u1, vec.Length);
                norm2 = NormVector(u1, RS);
            } while (Math.Abs(norm2 - norm1) > epsilon);
            vec = null;
            return norm2;
        }

        private void Eigensrch(double[] C, int L, double[] S, double[] U)
        {
            double[] T;
            double[] Q;

            double[] s;
            double[] a;
            double[] b;
            
            T = new double[L * L];
            Q = new double[L * CN];
            s = new double[CN * CN];

            a = new double[CN];
            b = new double[CN+1];

            MatrixTAndMatrixQIsResultOfLanczosTridiagonalizationMatrixA(C, L, L, T, Q);

            for (int i = 0; i < CN; i++)
            {
                a[i] = T[i * CN + i];
                b[i] = (i == 0) ? 0 : T[(i - 1) * CN + i];
            }
            b[CN] = 0;
/*
            for (int i = 0; i < CN; i++)
            {
                for (int j = 0; j < CN; j++)
                {
                    Console.Write("{0:f6}\t", T[i * CN + j]);
                }
                Console.Write("\n");
            }
            Console.Write("\n");
*/
            for (int t = 0; t < CN; t++)
            //Parallel.For(0, CN, t =>
            {
                //Gerschorins' theorem
                int evn = CN - t;
                double epsilon = 0.000001;

                double xmax = 0;
                double xmin = 0;

                for (int i = 0; i < CN; i++)
                {
                    double h = a[i] + Math.Abs(b[i]) + Math.Abs(b[i + 1]);
                    if (h > xmax) xmax = h;
                    h = a[i] - Math.Abs(b[i]) - Math.Abs(b[i + 1]);
                    if (h < xmin) xmin = h;
                }

                //Console.Write("[{0:f6}, {1:f6}] \n", xmin, xmax);

                do
                {
                    int NumOfEigensBefore = 0;
                    double d = 1.0;
                    double tmpMax = xmax;

                    xmax = (xmax + xmin) / 2.0;

                    for (int i = 0; i < CN; i++)
                    {
                        d = (a[i] - xmax) - (b[i] * b[i] / d);
                        if (d < 0) NumOfEigensBefore++;
                    }

                    if (NumOfEigensBefore < CN - t)
                    {
                        xmin = xmax;
                        xmax = tmpMax;
                    }

                    //Console.Write("xmax={0:f6} \n", xmax);
                    //Console.Write("xmin={0:f6} \n", xmin);
                } while (Math.Abs(xmax - xmin) > epsilon);
                S[t] = Math.Abs(xmax);

                double[] Tm;
                Tm = new double[CN * CN];
                for (int i = 0; i < CN; i++)
                    for (int j = 0; j < CN; j++)
                    {
                        if (i == j) Tm[i * CN + j] = T[i * CN + j] - S[t];
                        else
                            Tm[i * CN + j] = T[i * CN + j];
                    }

                //Eigen vector srch
                double[] ap;
                double[] bp;
                double[] c_coef;
                double[] qp;
                double[] cp;
                double[] dp;

                ap = new double[CN];
                bp = new double[CN];
                c_coef = new double[CN];

                qp = new double[CN * CN];
                cp = new double[CN * CN];
                dp = new double[CN * CN];

                //Invers Matrix Tm;
                GetApBp(Tm, ap, bp, CN);
                GetC_Сoef(c_coef, ap, bp, CN);
                GetQp(qp, ap, bp, c_coef, CN);
                GetCp(cp, ap, bp, c_coef, CN);
                GetDp(dp, ap, bp, c_coef, CN);

                //[cp]x[qp]x[dp]
                Array.Clear(Tm, 0, Tm.Length);
                double[] tmp_m;
                tmp_m = new double[CN * CN];
                for (int i = 0; i < CN; i++)
                {
                    for (int j = 0; j < CN; j++)
                    {
                        for (int k = 0; k < CN; k++)
                        {
                            tmp_m[i * CN + j] = tmp_m[i * CN + j] + cp[i * CN + k] * qp[k * CN + j];
                        }
                    }
                }

                for (int i = 0; i < CN; i++)
                {
                    for (int j = 0; j < CN; j++)
                    {
                        for (int k = 0; k < CN; k++)
                        {
                            Tm[i * CN + j] = Tm[i * CN + j] + tmp_m[i * CN + k] * dp[k * CN + j];
                        }
                    }
                }
/*
                for (int i = 0; i < CN; i++)
                {
                    for (int j = 0; j < CN; j++)
                    {
                        Console.Write("{0}\t", cp[i * CN + j]);
                    }
                    Console.Write("\n");
                }
                Console.Write("\n");
*/


                double[] us = new double[CN];
                for (int i = 0; i < CN; i++) us[i] = 1.0;
                double lambda;
                lambda = ComputeEigenU1Lambda1(Tm, CN, CN, us, epsilon);


                for (int i = 0; i < CN; i++)
                {
                    Console.Write("{0}\n", us[i]);
                }
                Console.Write("\n");
                Console.Write("\n");

                //});

            }


            for (int j = 0; j < CN; j++)
            {
                Console.Write("{0:f6}\n", S[j]);
            }
            Console.Write("\n");


            //TheQRfactorizationOfSymmetricTridiagonalMatrix(a, b, s, CN);
            //MatrixCThisResultOfMultiplyingMatrixAOnMatrixB(U, Q, L, CN, s, CN, CN);

            b = null;
        }
    }
}
