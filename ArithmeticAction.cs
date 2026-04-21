using System;
using System.Collections.Generic;
using System.Text;

namespace PolynomialExpressionResolver
{
    internal class ArithmeticAction
    {
        public ArithmeticActionType Type { get; set; }

        public double Perform(double x, double y)
        {
            switch (Type)
            {
                case ArithmeticActionType.Add:
                    return x + y;
                case ArithmeticActionType.Sub:
                    return x - y;
                case ArithmeticActionType.Mul:
                    return x * y;
                case ArithmeticActionType.Div:
                    return x / y;
                case ArithmeticActionType.Pow:
                    return Math.Pow(x, y);
            }

            throw new Exception($"{Type} is not supported");
        }
    }
}
