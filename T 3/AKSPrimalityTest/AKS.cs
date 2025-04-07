using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq; 

namespace AKSPrimalityTest
{
    public class Polynomial
    {
        public Dictionary<int, BigInteger> Coefficients { get; private set; }

        public Polynomial()
        {
            Coefficients = new Dictionary<int, BigInteger>();
        }

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

        public Polynomial Mod(int r, BigInteger n)
        {
            var resultCoeffs = new Dictionary<int, BigInteger>();
            foreach (var term in Coefficients)
            {
                int degree = term.Key;
                BigInteger coefficient = term.Value;

                BigInteger coeffModN = BigInteger.Remainder(coefficient, n);
                 if (coeffModN.Sign < 0) coeffModN += n; 

                if (coeffModN != BigInteger.Zero)
                {
                    int reducedDegree = degree % r;
                    resultCoeffs.TryGetValue(reducedDegree, out BigInteger existingCoeff);
                    BigInteger newCoeff = BigInteger.Remainder(existingCoeff + coeffModN, n);
                    if (newCoeff.Sign < 0) newCoeff += n;

                    if (newCoeff != BigInteger.Zero)
                    {
                        resultCoeffs[reducedDegree] = newCoeff;
                    }
                    else if (resultCoeffs.ContainsKey(reducedDegree))
                    {
                       resultCoeffs.Remove(reducedDegree); 
                    }
                }
            }
            return new Polynomial(resultCoeffs);
        }


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


                    int reducedDegree = resultDegreeRaw % r;
                    BigInteger coeffModN = BigInteger.Remainder(resultCoeffRaw, n);
                    if (coeffModN.Sign < 0) coeffModN += n; 


                    if (coeffModN != BigInteger.Zero)
                    {
                       resultCoeffs.TryGetValue(reducedDegree, out BigInteger existingCoeff);
                       BigInteger newCoeff = BigInteger.Remainder(existingCoeff + coeffModN, n);
                       if (newCoeff.Sign < 0) newCoeff += n; 

                        if (newCoeff != BigInteger.Zero)
                        {
                            resultCoeffs[reducedDegree] = newCoeff;
                        }
                         else if (resultCoeffs.ContainsKey(reducedDegree))
                        {
                            resultCoeffs.Remove(reducedDegree); 
                        }
                    }
                }
            }
            return new Polynomial(resultCoeffs);
        }

        public static Polynomial PowerMod(Polynomial poly, BigInteger exponent, int r, BigInteger n)
        {
            Polynomial result = new Polynomial(0, BigInteger.One);
            Polynomial basePoly = poly.Mod(r, n); 

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

        public bool EqualsMod(Polynomial other, int r, BigInteger n)
        {
            Polynomial thisMod = this.Mod(r, n);
            Polynomial otherMod = other.Mod(r, n);


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

             if (BigInteger.Pow(2, bits-1) == n && n > 1) return bits -1;

        }

        private static BigInteger BigIntSqrt(BigInteger n)
        {
            if (n < 0) throw new ArgumentOutOfRangeException(nameof(n), "Cannot compute square root of negative number.");
            if (n == 0) return 0;
            if (n < 4) return 1;

            BigInteger root = n >> (BigIntLog2(n)/2); 
            BigInteger lastRoot;
            do
            {
                lastRoot = root;
                root = (lastRoot + n / lastRoot) >> 1;
            } while (root < lastRoot);

            while (root * root > n) {
                 root--;
            }
            return root;
        }

         private static BigInteger NthRoot(BigInteger n, int b)
         {
             if (n < 0) throw new ArgumentException("N must be non-negative.");
             if (b <= 0) throw new ArgumentException("Exponent b must be positive.");
             if (n == 0) return 0;
             if (n == 1) return 1;
             if (b == 1) return n;
             if (b == 2) return BigIntSqrt(n);

             BigInteger low = 1;

             int logN = BigIntLog2(n);
             BigInteger high = BigInteger.One << (logN / b + 1); 

             BigInteger root = 0;

             while (low <= high)
             {
                 BigInteger mid = low + (high - low) / 2;
                 if (mid == 0) break; 

                 try {
                      BigInteger power = BigInteger.Pow(mid, b);
                      if (power == n) return mid;
                      if (power < n) {
                          root = mid; 
                          low = mid + 1;
                      } else {
                          high = mid - 1; 
                      }
                 } catch (OverflowException) {
                     high = mid - 1;
                 }
             }
             return root; 
         }

        public static bool IsPrime(BigInteger n)
        {
             // Handle trivial cases
             if (n <= 1) return false;
             if (n <= 3) return true; 
             if (n % 2 == 0 || n % 3 == 0) return false;


            if (IsPerfectPower(n))
            {

                return false; 
            }


            int log2n = BigIntLog2(n);
            int log2nSquared = log2n * log2n;


            BigInteger r = 2;
            while (true) 
            {
            
                if (BigInteger.GreatestCommonDivisor(r, n) > 1)
                {

                    return false; 
                }


                int order = ComputeOrder(r, n);

                if (order > log2nSquared)
                {

                    break; 
                }
                 if (r > n) 
                 {

                     return false; 
                 }
                r++;
            }
            
             for (BigInteger a_check = 2; a_check <= r; ++a_check) {
                 BigInteger gcd_a = BigInteger.GreatestCommonDivisor(a_check, n);
                 if (gcd_a > 1 && gcd_a < n) {

                     return false;
                 }
             }

            if (n <= r)
            {

                return true;
            }

             BigInteger phiR = EulerTotient(r);
             BigInteger sqrtPhiR = BigIntSqrt(phiR);
             BigInteger limit = sqrtPhiR * log2n; 




            for (BigInteger a = 1; a <= limit; a++)
            {
                var polyXplusA = new Polynomial(new Dictionary<int, BigInteger> { { 1, 1 }, { 0, a } });
                Polynomial leftPoly = Polynomial.PowerMod(polyXplusA, n, (int)r, n);

                 var polyXnPlusA_coeffs = new Dictionary<int, BigInteger> { { (int)(n % r), 1 }, { 0, a } }; 

                 if (n % r == 0) {
                      polyXnPlusA_coeffs[0] = BigInteger.Remainder(polyXnPlusA_coeffs[0] + 1, n);
                      if (polyXnPlusA_coeffs[0].Sign < 0) polyXnPlusA_coeffs[0] += n;
                      polyXnPlusA_coeffs.Remove((int)(n%r)); 
                       if (polyXnPlusA_coeffs[0] == 0 && polyXnPlusA_coeffs.ContainsKey(0)) polyXnPlusA_coeffs.Remove(0);
                 } else {

                      polyXnPlusA_coeffs[(int)(n % r)] = BigInteger.Remainder(polyXnPlusA_coeffs[(int)(n % r)], n);
                      if (polyXnPlusA_coeffs[(int)(n % r)].Sign < 0) polyXnPlusA_coeffs[(int)(n % r)] += n;
                      if (polyXnPlusA_coeffs[(int)(n % r)] == 0) polyXnPlusA_coeffs.Remove((int)(n % r));
                 }

                 if (polyXnPlusA_coeffs.ContainsKey(0)) {
                     polyXnPlusA_coeffs[0] = BigInteger.Remainder(polyXnPlusA_coeffs[0], n);
                     if (polyXnPlusA_coeffs[0].Sign < 0) polyXnPlusA_coeffs[0] += n;
                      if (polyXnPlusA_coeffs[0] == 0) polyXnPlusA_coeffs.Remove(0);
                 }


                Polynomial rightPoly = new Polynomial(polyXnPlusA_coeffs);



                if (!leftPoly.EqualsMod(rightPoly, (int)r, n))
                {

                    return false; 
                }
            }


            return true;
        }

        private static bool IsPerfectPower(BigInteger n)
        {
             if (n <= 3) return false; 


            int maxB = BigIntLog2(n);


            for (int b = 2; b <= maxB; b++)
            {

                BigInteger a = NthRoot(n, b);




                 if (a <= 1) {

                     break;
                 }



                try {
                     BigInteger checkPower = BigInteger.Pow(a, b);

                    if (checkPower == n)
                    {

                        return true;
                    }
                } catch (OverflowException) {


                    continue; 
                }

            }

            return false;
        }

        private static int ComputeOrder(BigInteger r_bi, BigInteger n)
        {
            if (r_bi > int.MaxValue) throw new ArgumentOutOfRangeException(nameof(r_bi), "Modulus r is too large for this implementation's order finding.");
            int r = (int)r_bi;
             if (r <= 1) return 0; 

            BigInteger nModR = n % r;
            if (nModR.Sign < 0) nModR += r;
             if (nModR == 0) return 0; 
             if (nModR == 1) return 1;

            BigInteger currentPower = nModR;
            int k = 1;

            while (k < r) 
            {
                k++;
                currentPower = (currentPower * nModR) % r;
                 if (currentPower.Sign < 0) currentPower += r;


                if (currentPower == 1)
                {
                    return k;
                }
            }
            return 0; 
        }

        private static BigInteger EulerTotient(BigInteger r)
        {
            if (r <= 0) return 0;
            BigInteger result = r;
            BigInteger p = 2;
            BigInteger tempR = r; 


            while (p * p <= tempR)
            {
                if (tempR % p == 0)
                {

                    result -= result / p;

                    while (tempR % p == 0)
                    {
                        tempR /= p;
                    }
                }
                p++;
            }


            if (tempR > 1)
            {
                result -= result / tempR;
            }

            return result;
        }
    }
}