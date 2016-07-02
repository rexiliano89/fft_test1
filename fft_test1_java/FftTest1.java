public class FftTest1 {
    public static void main(String args[]) {
        final int n = 8192; // data length (int), n >= 2, n = power of 2
        FFTsg fft = new FFTsg(n);
        
        double a[] = new double[n];
        
        for (int i = 0; i < n; i++) {
        	a[i] = Math.sin(2.0 * Math.PI * i / n * 10.0) + Math.cos(2.0 * Math.PI * i / n * 3.0);
        }
        
        double a_tmp0[] = a.clone();
        fft.rdft(1, a_tmp0);

        // M回測定
        final int M = 10000;
        
        long start = System.currentTimeMillis();

        for (int i = 0; i < M; i++) {
        	double a_tmp[] = a.clone();
        	fft.rdft(1, a_tmp);
        }
        
        long stop = System.currentTimeMillis();
        
        System.out.printf("time = %d ms\n", stop - start);
    }
}
