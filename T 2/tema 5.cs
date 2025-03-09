using System;

public class EulerTotientCalculator
{
    public static long EulerTotient(long n)
    {
        // Handle edge cases
        if (n <= 0) return 0;
        if (n == 1) return 1;

        long result = n;

        for (long p = 2; p * p <= n; p++)
        {
            if (n % p == 0)
            {
                while (n % p == 0)
                {
                    n /= p;
                }
                
                result -= result / p;
            }
        }


        if (n > 1)
        {
            result -= result / n;
        }

        return result;
    }
    public static bool VerifyEulerTotient()
    {
        bool test1 = EulerTotient(10) == 4;
        bool test2 = EulerTotient(36) == 12;

        bool test3 = EulerTotient(1) == 1;
        bool test4 = EulerTotient(0) == 0;

        bool test5 = EulerTotient(13) == 12;
        bool test6 = EulerTotient(123456) > 0;

        bool test7 = EulerTotient(17) == 16;

        return test1 && test2 && test3 && test4 && test5 && test6 && test7;
    }
    public static void Main()
    {
        // Verify the implementation
        Console.WriteLine($"Implementation Verified: {VerifyEulerTotient()}");

        // Demonstrate some example calculations
        long[] testNumbers = { 10, 36, 123456, 1000000007 };

        foreach (long num in testNumbers)
        {
            Console.WriteLine($"φ({num}) = {EulerTotient(num)}");
        }
    }
}