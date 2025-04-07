using System;
using System.Collections.Generic;

namespace AKSPrimalityTest
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("AKS Primality Test");
            Console.WriteLine("------------------");
            int[] testNumbers = { 2, 3, 4, 5, 7, 11, 13, 15, 17, 19, 23, 25, 29, 31,  1000003, 10000019,1000033, 1000037, 1000039, 1000041, 1000043, 1000049, 1000053, 1000061, 1000063 };
            Console.WriteLine("Testing the following numbers:");
            foreach (int n in testNumbers)
            {
                Console.Write(n + " ");
            }
            Console.WriteLine("\n\nResults:");
            foreach (int n in testNumbers)
            {
                Console.WriteLine($"Testing {n}...");
                string result = AKS.IsPrime(n) ? "prime" : "not prime";
                Console.WriteLine($"{n} is {result}");
            }

             // Initialize coefficients for AKS test

            Console.WriteLine("\nEnter a number to test (or 0 to exit):");
            while (true)
            {
                string input = Console.ReadLine();
                if (int.TryParse(input, out int number))
                {
                    if (number == 0)
                        break;

                    string result = AKS.IsPrime(number)? "prime" : "not prime";
                    Console.WriteLine($"{number} is {result}");
                    Console.WriteLine("\nEnter another number (or 0 to exit):");
                }
                else
                {
                    Console.WriteLine("Invalid input. Please enter a valid integer.");
                }
            }
        }
    }
}