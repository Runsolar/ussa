namespace pssaclass
{
    class pssaclass2
    {
        int N;
        int L;
        int K;

        double[] A;
        double[] R;
        double[] T;
        double[] Q;
        double[] y;
        double[] S;

        double[] AT;
        double[] X;

        public int ComponentsNum;

        public pssaclass2(double[] Signal, int Len, int Window, double[] Spectrum)
        {
            N = Len;
            L = Window;
            K = (N - L) + 1;

            A = new double[L * K];
            R = new double[L * K];
            T = new double[L * L];
            Q = new double[L * L];
            y = new double[L * L];
            S = new double[L];

            AT = new double[K * L];
            X = new double[L * L];

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
            //Step 3. T=Q*XQ
            this.MatrixTAndMatrixQIsResultOfLanczosTridiagonalizationMatrixA(X, L, L, T, Q);

            //Предварительный подсчет полного веса (ГК)
            double Ssum = 0;
            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < L; j++)
                {
                    Ssum = Ssum + T[i * L + j] * T[i * L + j];
                }
            }

            int MaxGnum = 0;
            double epsilon = 0.000001;
            double lambda;
            double[] Affect = new double[L];
            double[] u = new double[L];
            double[] v = new double[K];

            //Вычисление ГК
            do
            {
                Array.Clear(X, 0, X.Length);
                Array.Clear(R, 0, X.Length);

                for (int i = 0; i < L; i++) u[i] = 1.0;
                lambda = ComputeEigenU1Lambda1(T, L, L, u, epsilon);
                this.MulMatrixOnVector(T, L, L, u, v);

                for (int i = 0; i < L; i++)
                {
                    v[i] = v[i] / Math.Sqrt(lambda);
                }

                MatrixXThisResultOfMultiplyingVectorAOnVectorB(X, u, L, v, L);

                for (int i = 0; i < L; i++)
                {
                    for (int j = 0; j < L; j++)
                    {
                        X[(i * L + j)] = Math.Sqrt(lambda) * X[(i * L + j)];
                    }
                }

                this.MatrixAsubMatrixB(T, X, L, L);
                this.MulMatrixOnVector(Q, L, L, u, v);
                this.MulMatrixOnVector(AT, K, L, y, v);

                this.MatrixXThisResultOfMultiplyingVectorAOnVectorB(R, y, L, v, K);

                                for (int j = 0; j < L; j++)
                                {                             
                                    for (int i = 0; i < (j + 1); i++)
                                    {
                                        Spectrum[N * MaxGnum + j] = Spectrum[N * MaxGnum + j] + R[i * K + (j - i)];
                                    }
                                    Spectrum[N * MaxGnum + j] = Spectrum[N * MaxGnum + j] / (j + 1);

                                }

                                for (int j = L; j < K; j++)
                                {
                                    for (int i = 0; i < L; i++)
                                    {
                                        Spectrum[N * MaxGnum + j] = Spectrum[N * MaxGnum + j] + R[i * K + (j - i)];
                                    }
                                    Spectrum[N * MaxGnum + j] = Spectrum[N * MaxGnum + j] / L;
                                }

                                for (int j = K; j < N; j++)
                                {
                                    for (int i = 0; i < L - 1 - (j - K); i++)
                                    {
                                        Spectrum[N * MaxGnum + j] = Spectrum[N * MaxGnum + j] + R[(j - K + i) * K + (K - 1 - i)];
                                    }
                                    Spectrum[N * MaxGnum + j] = Spectrum[N * MaxGnum + j] / (N - j);
                                }

                Affect[MaxGnum] = (lambda * 100) / Ssum;
                MaxGnum++;
                //}while(Affect[MaxGnum-1] > 0.01);
 
            } while (MaxGnum <32);
            ComponentsNum = MaxGnum;
        }

        //Private method's    
        //X=(Ui)x(VTi)
        private void MatrixXThisResultOfMultiplyingVectorAOnVectorB(double[] X, double[] A, int L, double[] B, int K)
        {
            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    X[i * K + j] = A[i] * B[j];
                }
            }
        }

        private void MatrixAsubMatrixB(double[] A, double[] B, int RS, int CS)
        {
            for (int i = 0; i < RS; i++)
            {
                for (int j = 0; j < CS; j++)
                {
                    A[i * CS + j] = A[i * CS + j] - B[i * CS + j];
                }
            }
        }

        private void MulMatrixOnVector(double[] A, int RS, int CS, double[] u1, double[] vec)
        {
            //memset(vec, 0, RS*sizeof(double));
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
                //memcpy(u1, vec, RS*sizeof(double));
                //System.arraycopy(vec, 0, u1, 0, vec.length);

                Array.Copy(vec, u1, vec.Length);

                norm2 = NormVector(u1, RS);
            } while (Math.Abs(norm2 - norm1) > epsilon);
            vec = null;
            return norm2;
        }

        private void TransposeMatrixThisResultOfTransposeMatrix(double[] Matrix, int L, int K, double[] TransposeMatrix)
        {
            for (int i = 0; i < L; i++)
                for (int j = 0; j < K; j++)
                {
                    TransposeMatrix[j * L + i] = Matrix[i * K + j];
                }
        }

        private void MatrixCThisResultOfMultiplyingMatrixAOnMatrixB(double[] C, double[] A, int AN, int AM, double[] B, int BN, int BM)
        {
            //Arrays.fill(C, (double)0);
            Array.Clear(C, 0, C.Length);

            for (int i = 0; i < AN; i++)
            {
                for (int j = 0; j < BM; j++)
                {
                    for (int k = 0; k < BN; k++)
                    {
                        C[i * BM + j] = C[i * BM + j] + A[i * AM + k] * B[k * BM + j];
                    }
                }
            }
        }

        private void MatrixTAndMatrixQIsResultOfLanczosTridiagonalizationMatrixA(double[] A, int RS, int CS, double[] T, double[] Q)
        {
            double[] z = new double[RS];
            int m = 0;
            //memset(z, 0, RS*sizeof(double));
            //memset(T, 0, RS*CS*sizeof(double));
            //memset(Q, 0, RS*CS*sizeof(double));
            //Arrays.fill(z, (double)0);
            //Arrays.fill(T, (double)0);
            //Arrays.fill(Q, (double)0);
            Array.Clear(z, 0, z.Length);
            Array.Clear(T, 0, T.Length);
            Array.Clear(Q, 0, Q.Length);

            Q[0] = 1.0;
            for (int i = 0; i < CS; i++)
            {
                //VectorZisResultOfMultiplyingMatrixAOnVectorX(A, RS, CS, Q+i, z);
                for (int n = 0; n < RS; n++)
                {
                    z[n] = (double)0.0;
                    for (int j = 0; j < CS; j++)
                    {
                        z[n] = z[n] + A[n * CS + j] * Q[j * CS + i];
                    }
                }
                //
                //ArgumentCThisResultOfMultiplyingVectorZTonVectorX(T+CS*i+i, Q+i, z, RS);
                T[CS * i + i] = (double)0.0;
                for (int n = 0; n < RS; n++)
                {
                    T[CS * i + i] = T[CS * i + i] + Q[n * CS + i] * z[n];
                }
                //           
                if (i != CS - 1)
                {
                    for (int n = 0; n < RS; n++)
                    {
                        Q[n * CS + i + 1] = (i == 0) ? z[n] - T[CS * i + i] * Q[n * CS + i] : z[n] - T[CS * i + i] * Q[n * CS + i] - T[(i - 1) * CS + i] * Q[n * CS + i - 1];
                    }
                    double R;
                    for (int n = 0; n < m; n++)
                    {
                        R = 0;
                        for (int k = 0; k < RS; k++)
                        {
                            R = R + Q[k * CS + n] * Q[k * CS + i + 1];
                        }
                        for (int k = 0; k < RS; k++)
                        {
                            Q[k * CS + i + 1] = Q[k * CS + i + 1] - Q[k * CS + n] * R;
                        }
                    }
                    //ArgumentNormIsResultOfNormalizationVectorX(Q+i+1, CS, T+i*CS+i+1);
                    T[i * CS + i + 1] = 0.0;
                    for (int n = 0; n < RS; n++)
                        T[i * CS + i + 1] = T[i * CS + i + 1] + Q[n * CS + i + 1] * Q[n * CS + i + 1];
                    T[i * CS + i + 1] = Math.Sqrt(T[i * CS + i + 1]);
                    //
                    for (int n = 0; n < RS; n++)
                        Q[n * CS + i + 1] = Q[n * CS + i + 1] / T[i * CS + i + 1];
                    T[(i + 1) * CS + i] = T[i * CS + i + 1];
                }
                m++;
            }
            z = null;
        }

    }
}
