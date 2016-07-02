using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace fft_test1_cs
{
    class Program
    {
        const int n = 8192; // data length (int), n >= 2, n = power of 2

        static void Main(string[] args)
        {
            FFT4g fft = new FFT4g(n);

            double[] a = new double[n];

            for (int i = 0; i < n; i++)
            {
                a[i] = Math.Sin(2.0 * Math.PI * i / n * 10.0) + Math.Cos(2.0 * Math.PI * i / n * 3.0);
            }

            double[] a_tmp = new double[n];
            a.CopyTo(a_tmp, 0);
            fft.rdft(1, a_tmp);

            // M回測定
            const int M = 10000;

            System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            for (int i = 0; i < M; i++)
            {
                a.CopyTo(a_tmp, 0);
                fft.rdft(1, a_tmp);
            }

            sw.Stop();

            Console.WriteLine("time =" + sw.ElapsedMilliseconds + "ms");
        }
    }
}
