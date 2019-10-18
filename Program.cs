using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace pssaclass
{
    class Program
    {
        static void Main(string[] args)
        {
            int FRAME = 256;
            int ComponentsNum = 32;

            int Len = FRAME;
            int Window = FRAME / 2;
            //int Len = 32;
            //int Window = 4;

            int m = 0;

            double[] Vector = new double[Len];
            double[] Spectrum = new double[ComponentsNum * Len];

            string[] lines = System.IO.File.ReadAllLines(@args[0]);

            foreach (string line in lines)
            {
                if (m < Len)
                {
                    Vector[m] = Double.Parse(line);
                }
                m++;
            }

            //Singenerator
            double[] Sinus = new double[125*Len];
            int p=0;
            for(int i=0; i<125; i++)
            {
                for (int j = 0; j < Len; j++)
                    Sinus[i*Len+j]=Math.Sin(2 * Math.PI * (p + 100) * j * ((double)1 / 8192));
                                p=p+2;
            }

            p = 0;
            //double MFPE = 0;
            //for (int t = 0; t < 125; t++)
            {
               
                //Console.Write("Sinus {0} MGz\n", p + 100);
                
                /*
                            for (int i = 0; i < 125; i++)
                            {
                                for (int j = 0; j < 512; j++)
                                    Console.Write("{0:f6}\t", Sinus[i * Len + j]);
                                Console.WriteLine("");
                            }
                */

                //Random rnd = new Random();


                //for (int j = 0; j < Len; j++)
                //    Vector[j] = ((double)rnd.Next(-31620, 31620)) / 10000 + Sinus[t * Len + j];


                Stopwatch stopwatch = Stopwatch.StartNew();
                //pssaclass pssa = new pssaclass(Vector, Len, Window, Spectrum);
                //ExtremeSSG pssa = new ExtremeSSG(Vector, Len, Window, Spectrum, ComponentsNum);
                //ssgclass pssa = new ssgclass(Vector, Len, Window, Spectrum);
                BoostSSG pssa = new BoostSSG(Vector, Len, Window, Spectrum, ComponentsNum);

                stopwatch.Stop();

                Console.WriteLine("Time to SSG msec: {0}", stopwatch.ElapsedMilliseconds);
            
  /*                          Console.WriteLine("Vector:");
                            for (int i = 0; i < Len; i++)
                            {
                                Console.Write("{0:f6}\t", Vector[i]);
                            }
                            Console.WriteLine("");

                            Console.WriteLine("Spectrum:\n");
                            for (int i = 0; i < ComponentsNum; i++)
                            {
                                for (int j = 0; j < Len; j++)
                                {
                                    Console.Write("{0:f6}\t", Spectrum[i * Len + j]);
                                    //Console.Write("{0}\t", (Spectrum[i * Len + j] < 0) ? -1 : 1);
                    
                                }
                                Console.WriteLine("");
                            }
                            Console.WriteLine("");
 */            
/*
                double[] energy = new double[ComponentsNum];
                double[] frec = new double[ComponentsNum];

                for (int i = 0; i < ComponentsNum; i++)
                {
                    int n1 = 0;
                    int n2 = 0;
                    int n3 = 0;
                    int n = 0;

                    bool first = false;
                    bool start = true;

                    int[] candidate = new int[Len];

                    int MaxNum = 0;

                    for (int j = 1; j < Len - 1; j++)
                    {
                        double x1 = Spectrum[i * Len + j];
                        double x2 = Spectrum[i * Len + j + 1];

                        double y1 = Spectrum[i * Len + j - 1];
                        double y2 = Spectrum[i * Len + j];
                        double y3 = Spectrum[i * Len + j + 1];

                        if (x1 < 0 && x2 > 0)
                        {
                            candidate[j] = 1;
                        }
                        else
                        {
                            candidate[j] = 0;
                        }

                        if (y1 > 0 && y2 > 0 && y3 > 0 && ((y2 - y1) > 0) && ((y3 - y2) < 0))
                        {
                            energy[i] = energy[i] + Spectrum[i * Len + j];
                            MaxNum++;
                        }
                    }
                    energy[i] = energy[i] / MaxNum;

                    for (int j = 0; j < Len; j++)
                    {
                        if (candidate[j] == 1 && start)
                        {
                            if (first) start = false;

                            first = true;
                            n1 = j;
                        }

                        if (candidate[j] == 1)
                        {
                            n++;
                            n2 = n3;
                            n3 = j;
                        }

                    }

                    //frec[i] = (double)1.0 / (((n2 - n1) / (n - 3)) * (1.0 / 8192));
                }

                Array.Sort(frec);

                //Console.Write(" F={0:f6}\n", frec[2]);

                for (int i = 0; i < ComponentsNum-1; i++)
                {
                    if (90.0 < frec[i] && frec[i] < 370.0)
                    {
                        Console.Write(" F0={0:f6}\n", frec[i]);
                    }
                }
                Console.WriteLine("");


                //MFPE = MFPE + Math.Abs((p + 100) - (frec[0] + frec[1]) / 2) / (p + 100);

                p = p + 2;
 */ 
            }

           // MFPE = MFPE / 125;
           // Console.Write(" MFPE={0:f6}\n", MFPE);

        }

    }
}
