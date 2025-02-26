using System;

class BaseConverter
{
    private const string Digits = "0123456789ABCDEFGHIJKLMNOP";

    public static int ToDecimal(string number, int baseFrom)
    {
        int decimalValue = 0;
        foreach (char digit in number)
        {
            decimalValue = decimalValue * baseFrom + Digits.IndexOf(digit);
        }
        return decimalValue;
    }

    public static string FromDecimal(int number, int baseTo)
    {
        if (number == 0)
            return "0";
        
        string result = "";
        while (number > 0)
        {
            result = Digits[number % baseTo] + result;
            number /= baseTo;
        }
        return result;
    }

    public static string ConvertBase(string number, int baseFrom, int baseTo)
    {
        int decimalValue = ToDecimal(number, baseFrom);
        return FromDecimal(decimalValue, baseTo);
    }

    public static void Main()
    {
        string number = "1A"; 
        int baseFrom = 20;    
        int baseTo = 10;       
        Console.WriteLine(ConvertBase(number, baseFrom, baseTo));
    }
}
