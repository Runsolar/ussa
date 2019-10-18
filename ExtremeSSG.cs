using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace pssaclass
{
    class ExtremeSSG
    {
        int N;
        int L;
        int K;

        double[] A;
        double[] AT;
        double[] X;

        double[] S;
        double[] U;
        double[] VT;

        int CN;

        public ExtremeSSG(double[] Signal, int Len, int Window, double[] Spectrum, int ComponentsNum)
        {
            N = Len;
            L = Window;
            K = N - L + 1;

            CN = ComponentsNum;

            A = new double[L * K];
            AT = new double[K * L];
            X = new double[L * L];

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
            //Step 1. Transpose A -> AT
            this.TransposeMatrixThisResultOfTransposeMatrix(A, L, K, AT);
            //Step 2. X=A*AT
            this.MatrixCThisResultOfMultiplyingMatrixAOnMatrixB(X, A, L, K, AT, K, L);


            this.Eigensrch(X, L, S, U);

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
                        for (int k = 0; k < RS; k++){
                            R = R + Q[k * CN + n] * Q[k * CN + i + 1];
                        }
                        for (int k = 0; k < RS; k++){
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

        private void Eigensrch(double[] X, int L, double[] a, double[] U)
        {

            double[] T;
            double[] Q;

            double[] s;
            double[] b;

            T = new double[L * CN];
            Q = new double[L * CN];
            s = new double[CN * CN];


            b = new double[CN];

            MatrixTAndMatrixQIsResultOfLanczosTridiagonalizationMatrixA(X, L, L, T, Q);

            Array.Clear(a, 0, a.Length);
            Array.Clear(b, 0, b.Length);

            for (int i = 0; i < CN; i++)
            {
                a[i] = T[i * CN + i];
                if (i != CN - 1) b[i] = T[i * CN + i + 1];
            }

            TheQRfactorizationOfSymmetricTridiagonalMatrix(a, b, s, CN);
            MatrixCThisResultOfMultiplyingMatrixAOnMatrixB(U, Q, L, CN, s, CN, CN);

            b = null;

        }


        private void TheQRfactorizationOfSymmetricTridiagonalMatrix(double[] a, double[] b, double[] Qt, int RS)
        {

            int M = RS - 1;
            int m = 0;

            double r, t;

            double epsilon = 0;

            double[] c = new double[RS];
            double[] s = new double[RS];

            double p1 = 0;
            double p2 = 0;

            Array.Clear(Qt, 0, Qt.Length);

            for (int i = 0; i < RS; i++)
            {
                Qt[i * RS + i] = 1;
            }


            double D = 0;
            double S = 0;

            int[] mem = new int[RS];
            int g = 0;

            while (m<M)
            {
                //Сдвиг Вилкинсона
                //D = (a[m - 1] - a[m]) / 2;
                //S = (D == 0) ? a[m] - Math.Abs(b[m - 1]) : (a[m] - Math.Pow(b[m - 1], 2)) / (D + Math.Sign(D) * Math.Sqrt(Math.Pow(D, 2) + Math.Pow(b[m - 1], 2)));
        
                //сдвиг
                //D = (a[m + 1] - a[m]) / 2;                   
                //S = (D == 0) ? a[m + 1] - Math.Abs(b[m]) : (a[m + 1] - Math.Pow(b[m], 2)) / (D + Math.Sign(D) * Math.Sqrt(Math.Pow(D, 2) + Math.Pow(b[m], 2)));

                D = (a[g + 1] - a[g]) / 2;                   
                S = (D == 0) ? a[g + 1] - Math.Abs(b[g]) : (a[g + 1] - Math.Pow(b[g], 2)) / (D + Math.Sign(D) * Math.Sqrt(Math.Pow(D, 2) + Math.Pow(b[g], 2)));

                //Console.Write("{0:f6}\n", D);

                t = b[0];

                Array.Clear(s, 0, s.Length);
                Array.Clear(c, 0, c.Length);


                for (int i = 0; i < M; i++) a[i] = a[i] - S;

                for (int k = 0; k < M; k++)
                {
                    r = Math.Sqrt(a[k] * a[k] + t * t);

                    c[k] = a[k] / r;
                    s[k] = t / r;
                    a[k] = r;
                    t = b[k];
                    b[k] = t * c[k] + a[k + 1] * s[k];
                    a[k + 1] = -t * s[k] + a[k + 1] * c[k];

                    if (k < M - 1)
                    {
                        t = b[k + 1];
                        b[k + 1] = t * c[k];
                    }
                }

                //Accumulate H matrix in Qt
                for (int k = 0; k < M; k++)
                {
                    for (int j = 0; j < RS; j++)
                    {
                        p1 = Qt[j * RS + k];
                        p2 = Qt[j * RS + k + 1];

                        Qt[j * RS + k] = p1 * c[k] + p2 * s[k];
                        Qt[j * RS + k + 1] = -p1 * s[k] + p2 * c[k];
                    }
                }

                //RQ factorization
                //epsilon = 0.0;
                for (int k = 0; k < M; k++)
                {
                    a[k] = a[k] * c[k] + b[k] * s[k];
                    b[k] = a[k + 1] * s[k];
                    a[k + 1] = a[k + 1] * c[k];
                    //epsilon = epsilon + Math.Abs(b[k]);
                }
                for (int i = 0; i < M; i++) a[i] = a[i] + S;
                
                if (Math.Abs(b[m]) < 0.000001)
                {
                    m++;
                }
                
                g=(g==M-1)?0:g+1;

                if (Math.Abs(b[g]) < 0.000001)
                {
                    mem[g] = 1;
                }

                mem[M] = 0;
                for (int i = 0; i < M; i++)
                {
                    mem[M] = mem[M] + mem[i];
                }

                if ((double)mem[M] / M > 0.85) m = M;


                //Console.Write("{0:f6}\n", (double)mem[M] / M);
/*
                for (int i = 0; i < RS; i++)
                {

                    Console.Write("{0}\t", mem[i]);
                }
                Console.WriteLine("");
*/

                for (int i = 0; i < RS; i++)
                {

                    Console.Write("{0:f6}\t", a[i]);
                }
                Console.WriteLine("");

            } 

            c = null;
            s = null;
        }
        


/*
        private void TheQRfactorizationOfSymmetricTridiagonalMatrix(double[] a, double[] b, double[] Qt, int RS)
        {
            int m = RS - 1;
            double r, t;

            double epsilon = 0;

            double[] c = new double[RS];
            double[] s = new double[RS];

            double p1 = 0;
            double p2 = 0;

            Array.Clear(Qt, 0, Qt.Length);

            for (int i = 0; i < RS; i++)
            {
                Qt[i * RS + i] = 1;
            }

            do
            {

                for (int i = 0; i < RS; i++)
                {
                    Console.Write("{0:f6}\t", b[i]);
                }
                Console.WriteLine("");


                t = b[0];

                Array.Clear(s, 0, s.Length);
                Array.Clear(c, 0, c.Length);

                for (int k = 0; k < m; k++)
                {
                    r = Math.Sqrt(a[k] * a[k] + t * t);

                    c[k] = a[k] / r;
                    s[k] = t / r;
                    a[k] = r;
                    t = b[k];
                    b[k] = t * c[k] + a[k + 1] * s[k];
                    a[k + 1] = -t * s[k] + a[k + 1] * c[k];

                    if (k < m - 1)
                    {
                        t = b[k + 1];
                        b[k + 1] = t * c[k];
                    }
                }

                //Accumulate H matrix in Qt
                for (int k = 0; k < m; k++)
                {
                    for (int j = 0; j < RS; j++)
                    {
                        p1 = Qt[j * RS + k];
                        p2 = Qt[j * RS + k + 1];

                        Qt[j * RS + k] = p1 * c[k] + p2 * s[k];
                        Qt[j * RS + k + 1] = -p1 * s[k] + p2 * c[k];
                    }
                }
                //RQ factorization

                epsilon = 0.0;
                for (int k = 0; k < m; k++)
                //Parallel.For(0, m, k =>
                {
                    a[k] = a[k] * c[k] + b[k] * s[k];
                    b[k] = a[k + 1] * s[k];
                    a[k + 1] = a[k + 1] * c[k];
                    epsilon = epsilon + Math.Abs(b[k]);
                }

            } while (epsilon > 10.000000);

            c = null;
            s = null;
        }
*/

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
    }
}
