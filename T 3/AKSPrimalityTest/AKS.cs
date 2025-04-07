using System;
using System.Numerics;
using System.Collections.Generic; // For Polynomial Dictionary representation
using System.Linq; // For Polynomial operations

namespace AKSPrimalityTest
{
    /// <summary>
    /// Represents a polynomial with BigInteger coefficients.
    /// Uses a dictionary for sparse representation (degree -> coefficient).
    /// </summary>
    public class Polynomial
    {
        // Using a dictionary is efficient for sparse polynomials that occur
        // during modular exponentiation. Key = degree, Value = coefficient.
        public Dictionary<int, BigInteger> Coefficients { get; private set; }

        public Polynomial()
        {
            Coefficients = new Dictionary<int, BigInteger>();
        }

        // Constructor for simple polynomials like 'c' or 'X + a'
        public Polynomial(int degree, BigInteger coefficient) : this()
        {
            if (coefficient != BigInteger.Zero)
            {
                Coefficients[degree] = coefficient;
            }
        }

        public Polynomial(Dictionary<int, BigInteger> coeffs)
        {
             Coefficients = new Dictionary<int, BigInteger>(coeffs);
        }


        public BigInteger GetCoefficient(int degree)
        {
            return Coefficients.TryGetValue(degree, out BigInteger coeff) ? coeff : BigInteger.Zero;
        }

        public int Degree => Coefficients.Count == 0 ? -1 : Coefficients.Keys.Max();

        // Reduces the polynomial modulo (X^r - 1, n)
        // Effectively means:
        // 1. All coefficients are reduced modulo n.
        // 2. Any term X^k is replaced by X^(k mod r).
        public Polynomial Mod(int r, BigInteger n)
        {
            var resultCoeffs = new Dictionary<int, BigInteger>();
            foreach (var term in Coefficients)
            {
                int degree = term.Key;
                BigInteger coefficient = term.Value;

                BigInteger coeffModN = BigInteger.Remainder(coefficient, n);
                 if (coeffModN.Sign < 0) coeffModN += n; // Ensure positive remainder

                if (coeffModN != BigInteger.Zero)
                {
                    int reducedDegree = degree % r;
                    resultCoeffs.TryGetValue(reducedDegree, out BigInteger existingCoeff);
                    BigInteger newCoeff = BigInteger.Remainder(existingCoeff + coeffModN, n);
                    if (newCoeff.Sign < 0) newCoeff += n; // Ensure positive remainder

                    if (newCoeff != BigInteger.Zero)
                    {
                        resultCoeffs[reducedDegree] = newCoeff;
                    }
                    else if (resultCoeffs.ContainsKey(reducedDegree))
                    {
                       resultCoeffs.Remove(reducedDegree); // Remove zero terms
                    }
                }
            }
            return new Polynomial(resultCoeffs);
        }

        // Multiply two polynomials modulo (X^r - 1, n)
        public static Polynomial MultiplyMod(Polynomial p1, Polynomial p2, int r, BigInteger n)
        {
            var resultCoeffs = new Dictionary<int, BigInteger>();

            foreach (var term1 in p1.Coefficients)
            {
                foreach (var term2 in p2.Coefficients)
                {
                    int degree1 = term1.Key;
                    BigInteger coeff1 = term1.Value;
                    int degree2 = term2.Key;
                    BigInteger coeff2 = term2.Value;

                    int resultDegreeRaw = degree1 + degree2;
                    BigInteger resultCoeffRaw = coeff1 * coeff2;

                    // Reduce modulo (X^r - 1, n)
                    int reducedDegree = resultDegreeRaw % r;
                    BigInteger coeffModN = BigInteger.Remainder(resultCoeffRaw, n);
                    if (coeffModN.Sign < 0) coeffModN += n; // Ensure positive remainder


                    if (coeffModN != BigInteger.Zero)
                    {
                       resultCoeffs.TryGetValue(reducedDegree, out BigInteger existingCoeff);
                       BigInteger newCoeff = BigInteger.Remainder(existingCoeff + coeffModN, n);
                       if (newCoeff.Sign < 0) newCoeff += n; // Ensure positive remainder

                        if (newCoeff != BigInteger.Zero)
                        {
                            resultCoeffs[reducedDegree] = newCoeff;
                        }
                         else if (resultCoeffs.ContainsKey(reducedDegree))
                        {
                            resultCoeffs.Remove(reducedDegree); // Remove zero terms
                        }
                    }
                }
            }
            return new Polynomial(resultCoeffs);
        }

        // Computes poly^exponent mod (X^r - 1, n) using binary exponentiation
        public static Polynomial PowerMod(Polynomial poly, BigInteger exponent, int r, BigInteger n)
        {
            // Identity polynomial is 1 (degree 0, coefficient 1)
            Polynomial result = new Polynomial(0, BigInteger.One);
            Polynomial basePoly = poly.Mod(r, n); // Ensure base is reduced initially

            while (exponent > 0)
            {
                if (exponent % 2 == 1)
                {
                    result = MultiplyMod(result, basePoly, r, n);
                }
                basePoly = MultiplyMod(basePoly, basePoly, r, n);
                exponent /= 2;
            }
            return result;
        }

        // Checks if two polynomials are equal modulo (X^r - 1, n)
        public bool EqualsMod(Polynomial other, int r, BigInteger n)
        {
            // Reduce both polynomials first to ensure fair comparison
            Polynomial thisMod = this.Mod(r, n);
            Polynomial otherMod = other.Mod(r, n);

            // Check if they have the same terms (degree and coefficient)
            if (thisMod.Coefficients.Count != otherMod.Coefficients.Count)
                return false;

            foreach (var term in thisMod.Coefficients)
            {
                if (!otherMod.Coefficients.TryGetValue(term.Key, out BigInteger otherCoeff) || term.Value != otherCoeff)
                {
                    return false;
                }
            }
            return true;
        }

        public override string ToString()
        {
            if (Coefficients.Count == 0) return "0";
            var sortedTerms = Coefficients.OrderByDescending(kvp => kvp.Key);
            return string.Join(" + ", sortedTerms.Select(kvp => $"{kvp.Value}*X^{kvp.Key}"));
        }
    }


    public class AKS
    {
        // Helper function for BigInteger Logarithm (base 2)
        private static int BigIntLog2(BigInteger n)
        {
            if (n <= 0) throw new ArgumentOutOfRangeException(nameof(n), "Logarithm argument must be positive.");
            byte[] bytes = n.ToByteArray();
            int bits = (bytes.Length - 1) * 8;
            byte topByte = bytes[bytes.Length - 1];
            while (topByte > 0)
            {
                topByte >>= 1;
                bits++;
            }
            // Adjust if n is exactly a power of 2 (optional, usually log floor is fine)
             if (BigInteger.Pow(2, bits-1) == n && n > 1) return bits -1;
            return bits > 0 ? bits -1 : 0; // Effectively floor(log2(n))
        }

         // Helper function for integer square root using binary search
        private static BigInteger BigIntSqrt(BigInteger n)
        {
            if (n < 0) throw new ArgumentOutOfRangeException(nameof(n), "Cannot compute square root of negative number.");
            if (n == 0) return 0;
            if (n < 4) return 1;

            BigInteger root = n >> (BigIntLog2(n)/2); // Initial guess using bit shifts for approx sqrt
            BigInteger lastRoot;
            do
            {
                lastRoot = root;
                root = (lastRoot + n / lastRoot) >> 1; // Newton's method step (integer division)
            } while (root < lastRoot); // Iterate until convergence

            // It might slightly overshoot, check
            while (root * root > n) {
                 root--;
            }
            return root;
        }

        // Helper: Nth root using binary search
         private static BigInteger NthRoot(BigInteger n, int b)
         {
             if (n < 0) throw new ArgumentException("N must be non-negative.");
             if (b <= 0) throw new ArgumentException("Exponent b must be positive.");
             if (n == 0) return 0;
             if (n == 1) return 1;
             if (b == 1) return n;
             if (b == 2) return BigIntSqrt(n);

             BigInteger low = 1;
             // An upper bound can be estimated. If a^b = n, then a <= n. A tighter bound: a <= 2^(ceil(log2(n)/b))
             int logN = BigIntLog2(n);
             BigInteger high = BigInteger.One << (logN / b + 1); // 2^(logN/b + 1)

             BigInteger root = 0;

             while (low <= high)
             {
                 BigInteger mid = low + (high - low) / 2;
                 if (mid == 0) break; // Avoid issues with mid=0

                 try {
                      BigInteger power = BigInteger.Pow(mid, b);
                      if (power == n) return mid;
                      if (power < n) {
                          root = mid; // Possible candidate, try higher
                          low = mid + 1;
                      } else {
                          high = mid - 1; // Too high
                      }
                 } catch (OverflowException) {
                     // If BigInteger.Pow overflows, mid is definitely too high
                     high = mid - 1;
                 }
             }
             return root; // Returns floor(n^(1/b))
         }


        /// <summary>
        /// Main AKS primality test function using BigInteger
        /// </summary>
        public static bool IsPrime(BigInteger n)
        {
             // Handle trivial cases
             if (n <= 1) return false;
             if (n <= 3) return true; 
             if (n % 2 == 0 || n % 3 == 0) return false;


            if (IsPerfectPower(n))
            {

                return false; // Composites
            }


            int log2n = BigIntLog2(n);
            int log2nSquared = log2n * log2n;


            BigInteger r = 2;
            while (true) // We are guaranteed to find such an r
            {
                // Ensure r and n are coprime for order calculation
                if (BigInteger.GreatestCommonDivisor(r, n) > 1)
                {

                    return false; 
                }


                int order = ComputeOrder(r, n);

                if (order > log2nSquared)
                {

                    break; // Found our r
                }
                 if (r > n) // Safety break, though theoretically not needed if n > 2
                 {

                     return false; // Or handle as prime if IsPerfectPower was passed? AKS theorem edge case. Usually n<=r step handles this.
                 }
                r++;
            }
            
             for (BigInteger a_check = 2; a_check <= r; ++a_check) {
                 BigInteger gcd_a = BigInteger.GreatestCommonDivisor(a_check, n);
                 if (gcd_a > 1 && gcd_a < n) {

                     return false;
                 }
             }


            // --- Step 4: If n <= r, output prime ---
            // Note: Using BigInteger for comparison
            if (n <= r)
            {

                return true;
            }

             BigInteger phiR = EulerTotient(r);
             BigInteger sqrtPhiR = BigIntSqrt(phiR);
             BigInteger limit = sqrtPhiR * log2n; 




            for (BigInteger a = 1; a <= limit; a++)
            {
                // Construct LHS polynomial: (X + a)^n
                // Represent (X + a) as {1: 1, 0: a} -> coeff for X^1 is 1, coeff for X^0 is a
                var polyXplusA = new Polynomial(new Dictionary<int, BigInteger> { { 1, 1 }, { 0, a } });
                Polynomial leftPoly = Polynomial.PowerMod(polyXplusA, n, (int)r, n); // r is int, n is BigInt

                // Construct RHS polynomial: X^n + a
                // Represent X^n + a as {n: 1, 0: a}
                 // IMPORTANT: Reduce this modulo (X^r - 1, n) *before* comparison
                 var polyXnPlusA_coeffs = new Dictionary<int, BigInteger> { { (int)(n % r), 1 }, { 0, a } }; // X^n becomes X^(n mod r)
                 // Need to handle potential coefficient addition if n % r == 0
                 if (n % r == 0) {
                      polyXnPlusA_coeffs[0] = BigInteger.Remainder(polyXnPlusA_coeffs[0] + 1, n);
                      if (polyXnPlusA_coeffs[0].Sign < 0) polyXnPlusA_coeffs[0] += n;
                      polyXnPlusA_coeffs.Remove((int)(n%r)); // Remove the X^n term, it merged with constant
                       if (polyXnPlusA_coeffs[0] == 0 && polyXnPlusA_coeffs.ContainsKey(0)) polyXnPlusA_coeffs.Remove(0);
                 } else {
                     // Ensure coefficient for X^n (which is X^(n mod r)) is reduced mod n
                      polyXnPlusA_coeffs[(int)(n % r)] = BigInteger.Remainder(polyXnPlusA_coeffs[(int)(n % r)], n);
                      if (polyXnPlusA_coeffs[(int)(n % r)].Sign < 0) polyXnPlusA_coeffs[(int)(n % r)] += n;
                      if (polyXnPlusA_coeffs[(int)(n % r)] == 0) polyXnPlusA_coeffs.Remove((int)(n % r));
                 }
                 // Ensure constant term 'a' is reduced mod n
                 if (polyXnPlusA_coeffs.ContainsKey(0)) {
                     polyXnPlusA_coeffs[0] = BigInteger.Remainder(polyXnPlusA_coeffs[0], n);
                     if (polyXnPlusA_coeffs[0].Sign < 0) polyXnPlusA_coeffs[0] += n;
                      if (polyXnPlusA_coeffs[0] == 0) polyXnPlusA_coeffs.Remove(0);
                 }


                Polynomial rightPoly = new Polynomial(polyXnPlusA_coeffs);


                // Compare LHS and RHS modulo (X^r - 1, n)
                if (!leftPoly.EqualsMod(rightPoly, (int)r, n))
                {

                    return false; 
                }
            }


            return true;
        }

        private static bool IsPerfectPower(BigInteger n)
        {
             if (n <= 3) return false; // 1, 2, 3 are not perfect powers in this context

            // Check up to b = log2(n) because if n = a^b, then a >= 2, so 2^b <= a^b = n, implies b <= log2(n)
            int maxB = BigIntLog2(n);


            for (int b = 2; b <= maxB; b++)
            {
                // Find the integer b-th root of n
                BigInteger a = NthRoot(n, b);



                // If a <= 1, further checks with larger b won't yield a > 1
                 if (a <= 1) {

                     break;
                 }


                // Check if a^b actually equals n
                try {
                     BigInteger checkPower = BigInteger.Pow(a, b);

                    if (checkPower == n)
                    {

                        return true;
                    }
                } catch (OverflowException) {
                    // This shouldn't happen if NthRoot is correct, but as a safeguard

                    continue; // Or potentially break if 'a' seems reasonable?
                }

            }

            return false;
        }

        private static int ComputeOrder(BigInteger r_bi, BigInteger n)
        {
            // The modulus 'r' in ord_r(n) is typically small enough to be an int.
            // Let's assume r fits in int for the loop limit, but use BigInteger n.
            if (r_bi > int.MaxValue) throw new ArgumentOutOfRangeException(nameof(r_bi), "Modulus r is too large for this implementation's order finding.");
            int r = (int)r_bi;
             if (r <= 1) return 0; // Order is not defined for r <= 1

             // We assume GCD(r, n) == 1 was checked before calling this function.

            // Calculate n mod r first
            BigInteger nModR = n % r;
            if (nModR.Sign < 0) nModR += r;
             if (nModR == 0) return 0; // Should not happen if GCD is 1, unless r=1
             if (nModR == 1) return 1; // Order is 1

            BigInteger currentPower = nModR;
            int k = 1;

            // Maximum order is phi(r), which is <= r-1.
            // Loop up to r as a safe upper bound.
            while (k < r) // Or k <= phi(r) if calculated
            {
                k++;
                currentPower = (currentPower * nModR) % r;
                 if (currentPower.Sign < 0) currentPower += r;


                if (currentPower == 1)
                {
                    return k;
                }
            }
            return 0; // Indicate failure or non-existence under assumption GCD=1
        }

        private static BigInteger EulerTotient(BigInteger r)
        {
            if (r <= 0) return 0;
            BigInteger result = r;
            BigInteger p = 2;
            BigInteger tempR = r; // Work with a temporary copy

            // Iterate through potential prime factors up to sqrt(r)
            while (p * p <= tempR)
            {
                if (tempR % p == 0)
                {
                    // Subtract multiples of p
                    result -= result / p;
                    // Remove all factors of p from tempR
                    while (tempR % p == 0)
                    {
                        tempR /= p;
                    }
                }
                p++;
            }

            // If tempR > 1 at this point, it must be a prime factor itself
            if (tempR > 1)
            {
                result -= result / tempR;
            }

            return result;
        }

        // Note: The original code had GCD(int, int). We need GCD(BigInteger, BigInteger).
        // System.Numerics.BigInteger has BigInteger.GreatestCommonDivisor(a, b).

    }
}