using System.Numerics;

namespace AKSPrimalityTest
{
    public class AKS
    {
        /// <summary>
        /// Main AKS primality test function
        /// </summary>
        public static string IsPrime(int n)
        {
            // Step 1: Check if n is a perfect power
            if (IsPerfectPower(n))
            {
                return "composite";
            }

            // Step 2: Find the smallest r such that ord_r(n) > (log2(n))^2
            int r = FindSmallestR(n);
            if (r == 0) // If r and n are not coprime
            {
                return "composite";
            }

            // Step 3: Trial division by all a ≤ min(r, n-1)
            for (int a = 2; a < Math.Min(r, n); a++)
            {
                if (n % a == 0)
                {
                    return "composite";
                }
            }

            // Step 4: If n ≤ r, output prime
            if (n <= r)
            {
                return "prime";
            }

            // Step 5: Polynomial check
            if (PolynomialCheck(n, r))
            {
                return "prime";
            }
            else
            {
                return "composite";
            }
        }

        /// <summary>
        /// Check if n is a perfect power (n = a^b for integers a > 1 and b > 1)
        /// </summary>
        private static bool IsPerfectPower(int n)
        {
            int maxLog = (int)Math.Log2(n);
            for (int b = 2; b <= maxLog; b++)
            {
                double aDouble = Math.Pow(n, 1.0 / b);
                int a = (int)Math.Round(aDouble);

                // Check if a^b equals n
                if (BigInteger.Pow(a, b) == n)
                {
                    return true;
                }
            }
            return false;
        }

        /// <summary>
        /// Compute the multiplicative order of n modulo r
        /// (smallest positive integer k such that n^k ≡ 1 (mod r))
        /// </summary>
        private static int ComputeOrder(int r, int n)
        {
            // First check if r and n are coprime
            if (GCD(r, n) > 1)
            {
                return 0; // No multiplicative order exists if r and n are not coprime
            }

            int k = 1;
            long power = r % n; // Start with r mod n

            while (power != 1)
            {
                power = (power * r) % n; // Multiply by r each time
                k++;

                // If we've gone through all possible remainders without finding 1
                if (k > n)
                {
                    return 0; // This shouldn't happen if r and n are coprime
                }
            }

            return k;
        }

        /// <summary>
        /// Find the smallest r such that ord_r(n) > (log2(n))^2
        /// </summary>
        private static int FindSmallestR(int n)
        {
            double log2nSquared = Math.Pow(Math.Log2(n), 2);
            int r = 2;
            while (r < n)
            {
                int order = ComputeOrder(r, n);
                if (order > log2nSquared)
                {
                    return r;
                }
                r++;
            }
            return r;
        }

        /// <summary>
        /// Calculate Euler's totient function φ(n)
        /// </summary>
        private static int EulerTotient(int n)
        {
            int result = n;
            int p = 2;
            while (p * p <= n)
            {
                if (n % p == 0)
                {
                    while (n % p == 0)
                    {
                        n /= p;
                    }
                    result -= result / p;
                }
                p++;
            }
            if (n > 1)
            {
                result -= result / n;
            }
            return result;
        }

        /// <summary>
        /// Check if (X+a)^n ≡ X^n+a (mod X^r-1, n) for required values of a
        /// </summary>
        private static bool PolynomialCheck(int n, int r)
        {
            int phiR = EulerTotient(r);
            int limit = (int)(Math.Sqrt(phiR) * Math.Log2(n));

            for (int a = 1; a <= limit; a++)
            {
                if (GCD(a, n) > 1)
                {
                    return false; // n is composite if gcd(a, n) > 1
                }

                // Create polynomials for (X+a) and (X^n+a)
                int[] leftPoly = new int[] { a, 1 }; // X + a
                int[] rightPoly = new int[n + 1];
                rightPoly[0] = a; // constant term
                rightPoly[n] = 1; // X^n coefficient

                // Calculate (X+a)^n mod (X^r-1, n)
                int[] leftResult = PolyPowerMod(leftPoly, n, r, n);

                // Reduce right_poly modulo (X^r-1, n)
                int[] rightResult = PolyMod(rightPoly, r, n);

                if (!PolyEqual(leftResult, rightResult, r, n))
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Reduce polynomial modulo (X^r-1, modulus)
        /// </summary>
        private static int[] PolyMod(int[] poly, int r, int modulus)
        {
            int[] result = new int[r];

            // For each term in the polynomial
            for (int i = 0; i < poly.Length; i++)
            {
                // Compute i % r to fold higher degree terms back
                int idx = i % r;
                result[idx] = (result[idx] + poly[i]) % modulus;
            }

            return result;
        }

        /// <summary>
        /// Multiply two polynomials and reduce modulo (X^r-1, modulus)
        /// </summary>
        private static int[] PolyMultiply(int[] poly1, int[] poly2, int r, int modulus)
        {
            // Initialize the result array with zeros
            int[] result = new int[r];

            // For each term in the first polynomial
            for (int i = 0; i < poly1.Length; i++)
            {
                // For each term in the second polynomial
                for (int j = 0; j < poly2.Length; j++)
                {
                    // Compute the position in the result (folded by X^r-1)
                    int idx = (i + j) % r;
                    // Add the product of the coefficients
                    result[idx] = (result[idx] + (poly1[i] * poly2[j]) % modulus) % modulus;
                }
            }

            return result;
        }

        /// <summary>
        /// Compute polynomial^exponent mod (X^r-1, modulus) using binary exponentiation
        /// </summary>
        private static int[] PolyPowerMod(int[] poly, int exponent, int r, int modulus)
        {
            // Base polynomial = 1 (represented as a polynomial)
            int[] result = new int[r];
            result[0] = 1;

            // Ensure poly has length r
            int[] basePoly = PolyMod(poly, r, modulus);

            // Binary exponentiation
            while (exponent > 0)
            {
                if (exponent % 2 == 1)
                {
                    result = PolyMultiply(result, basePoly, r, modulus);
                }

                basePoly = PolyMultiply(basePoly, basePoly, r, modulus);
                exponent /= 2;
            }

            return result;
        }

        /// <summary>
        /// Check if two polynomials are equal in the ring R = (Z/nZ)[X]/(X^r-1)
        /// </summary>
        private static bool PolyEqual(int[] poly1, int[] poly2, int r, int modulus)
        {
            // Ensure both polynomials are reduced modulo (X^r-1, n)
            int[] p1 = PolyMod(poly1, r, modulus);
            int[] p2 = PolyMod(poly2, r, modulus);

            // Compare each coefficient
            for (int i = 0; i < r; i++)
            {
                if (p1[i] != p2[i])
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Calculate the greatest common divisor of two integers
        /// </summary>
        private static int GCD(int a, int b)
        {
            while (b != 0)
            {
                int temp = b;
                b = a % b;
                a = temp;
            }
            return a;
        }
    }
}