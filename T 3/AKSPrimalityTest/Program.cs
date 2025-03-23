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

            int[] testNumbers = { 2, 3, 4, 5, 7, 11, 13, 15, 17, 19, 23, 25, 29, 31, (int)Math.Pow(2, 24) + 1 };
            foreach (int n in testNumbers)
            {
                string result = AKS.IsPrime(n);
                Console.WriteLine($"{n} is {result}");
            }

            Console.WriteLine("\nEnter a number to test (or 0 to exit):");
            while (true)
            {
                string input = Console.ReadLine();
                if (int.TryParse(input, out int number))
                {
                    if (number == 0)
                        break;

                    string result = AKS.IsPrime(number);
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