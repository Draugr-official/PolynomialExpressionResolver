using System;
using System.Collections.Generic;
using System.Linq;

namespace PolynomialExpressionResolver
{
    internal class Program
    {

        static double[] xAxisLinear = { 0, 1, 2, 3, 4, 5, 6, 7 };
        static double[] yAxisLinear = { 20, 25, 30, 35, 40, 45, 50, 55 };

        static double[] xAxisPolynomial = { 13, 14, 15 };
        static double[] yAxisPolynomial = { 5.63, 6, 5.63 };


        static double[] xAxisPolynomial2ndGrade = { 0, 1, 2, 3, 4, 5, 6, 7 };
        static double[] yAxisPolynomial2ndGrade = { 110, 110, 96, 93, 88, 71, 70, 49 };

        static void Main(string[] args)
        {
            double[] xAxis = xAxisPolynomial;
            double[] yAxis = yAxisPolynomial;

            var poly = FindPolynomialExpression(xAxis, yAxis, tolerance: 1e-7);
            if (poly != null)
            {
                Console.WriteLine("Polynomial fit found:");
                Console.WriteLine(poly.Replace(",", "."));
            }
            else
            {
                var matches = FindExpressions(xAxis, yAxis, maxDepth: 3, minConst: -20, maxConst: 20, tolerance: 1e-7);

                if (matches.Count == 0)
                {
                    Console.WriteLine("No matching expression found.");
                }
                else
                {
                    Console.WriteLine("Found expressions:");
                    foreach (var m in matches.Take(20))
                        Console.WriteLine(m.Replace(",", "."));
                }
            }
        }

        static string? FindPolynomialExpression(double[] x, double[] y, double tolerance)
        {
            if (x.Length != y.Length || x.Length == 0) return null;
            int n = x.Length;

            if (n >= 3)
            {
                var ls = SolveLeastSquaresQuadratic(x, y);
                if (ls != null)
                    return FormatPolynomial(ls, tolerance);
            }

            for (int degree = 0; degree <= n - 1; degree++)
            {
                int m = degree + 1;


                double[]? coeffs = null;

                var idx = new int[m];
                bool found = false;

                void TryCombination(int pos, int start)
                {
                    if (found) return;
                    if (pos == m)
                    {
                        double[,] V = new double[m, m];
                        double[] b = new double[m];
                        for (int ii = 0; ii < m; ii++)
                        {
                            double xp = 1.0;
                            double xv = x[idx[ii]];
                            for (int j = 0; j < m; j++)
                            {
                                V[ii, j] = xp;
                                xp *= xv;
                            }
                            b[ii] = y[idx[ii]];
                        }

                        var sol = GaussianSolve(V, b);
                        if (sol != null)
                        {
                            coeffs = sol;
                            found = true;
                        }
                        return;
                    }

                    for (int i = start; i <= n - (m - pos); i++)
                    {
                        idx[pos] = i;
                        TryCombination(pos + 1, i + 1);
                        if (found) return;
                    }
                }

                TryCombination(0, 0);
                if (coeffs == null)
                {
                    double[,] normal = new double[m, m];
                    double[] rhs = new double[m];
                    for (int i = 0; i < m; i++)
                    {
                        for (int j = 0; j < m; j++)
                        {
                            double sum = 0.0;
                            for (int k = 0; k < n; k++)
                                sum += Math.Pow(x[k], i + j);
                            normal[i, j] = sum;
                        }

                        double s = 0.0;
                        for (int k = 0; k < n; k++)
                            s += y[k] * Math.Pow(x[k], i);
                        rhs[i] = s;
                    }

                    coeffs = GaussianSolve(normal, rhs);
                }
                if (coeffs == null) continue;

                bool ok = true;
                for (int k = 0; k < n; k++)
                {
                    double val = 0.0;
                    double xp = 1.0;
                    for (int j = 0; j < m; j++)
                    {
                        val += coeffs[j] * xp;
                        xp *= x[k];
                    }
                    if (double.IsNaN(val) || double.IsInfinity(val) || Math.Abs(val - y[k]) > tolerance)
                    {
                        ok = false;
                        break;
                    }
                }

                if (ok)
                    return FormatPolynomial(coeffs, tolerance);
            }

            return null;
        }

        static string FormatPolynomial(double[] coeffs, double tol)
        {
            var parts = new List<string>();
            for (int i = coeffs.Length - 1; i >= 0; i--)
            {
                double c = coeffs[i];
                if (Math.Abs(c) < tol) continue;

                string coeffStr;
                if (Math.Abs(c - Math.Round(c)) < tol)
                    coeffStr = Math.Round(c).ToString();
                else
                    coeffStr = c.ToString("G6");

                if (i == 0)
                {
                    parts.Add(coeffStr);
                }
                else if (i == 1)
                {
                    if (coeffStr == "1") parts.Add("x");
                    else if (coeffStr == "-1") parts.Add("-x");
                    else parts.Add($"{coeffStr}*x");
                }
                else
                {
                    if (coeffStr == "1") parts.Add($"x^{i}");
                    else if (coeffStr == "-1") parts.Add($"-x^{i}");
                    else parts.Add($"{coeffStr}*x^{i}");
                }
            }

            if (parts.Count == 0) return "0";

            string expr = parts[0];
            for (int i = 1; i < parts.Count; i++)
            {
                if (parts[i].StartsWith("-")) expr += " - " + parts[i].Substring(1);
                else expr += " + " + parts[i];
            }
            return expr;
        }

        static double[]? GaussianSolve(double[,] Aorig, double[] borig)
        {
            int n = borig.Length;
            double[,] A = new double[n, n];
            double[] b = new double[n];
            Array.Copy(borig, b, n);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    A[i, j] = Aorig[i, j];

            for (int i = 0; i < n; i++)
            {
                int pivot = i;
                double max = Math.Abs(A[i, i]);
                for (int r = i + 1; r < n; r++)
                {
                    double v = Math.Abs(A[r, i]);
                    if (v > max) { max = v; pivot = r; }
                }
                if (Math.Abs(A[pivot, i]) < 1e-12) return null;

                if (pivot != i)
                {
                    for (int c = i; c < n; c++) { var tmp = A[i, c]; A[i, c] = A[pivot, c]; A[pivot, c] = tmp; }
                    var tb = b[i]; b[i] = b[pivot]; b[pivot] = tb;
                }

                double diag = A[i, i];
                for (int c = i; c < n; c++) A[i, c] /= diag;
                b[i] /= diag;

                for (int r = 0; r < n; r++)
                {
                    if (r == i) continue;
                    double factor = A[r, i];
                    if (Math.Abs(factor) < 1e-12) continue;
                    for (int c = i; c < n; c++) A[r, c] -= factor * A[i, c];
                    b[r] -= factor * b[i];
                }
            }

            return b;
        }

        static double[]? SolveLeastSquaresQuadratic(double[] x, double[] y)
        {
            int n = x.Length;

            double[,] V = new double[n, 3];
            for (int i = 0; i < n; i++)
            {
                V[i, 0] = 1.0;
                V[i, 1] = x[i];
                V[i, 2] = x[i] * x[i];
            }

            double[,] Q = new double[n, 3];
            double[,] R = new double[3, 3];
            for (int k = 0; k < 3; k++)
            {
                double[] v = new double[n];
                for (int i = 0; i < n; i++) v[i] = V[i, k];

                for (int j = 0; j < k; j++)
                {
                    double r = 0.0;
                    for (int i = 0; i < n; i++) r += Q[i, j] * v[i];
                    R[j, k] = r;
                    for (int i = 0; i < n; i++) v[i] -= r * Q[i, j];
                }

                double norm = 0.0;
                for (int i = 0; i < n; i++) norm += v[i] * v[i];
                norm = Math.Sqrt(norm);
                if (norm < 1e-15) return null;
                R[k, k] = norm;
                for (int i = 0; i < n; i++) Q[i, k] = v[i] / norm;
            }

            double[] qty = new double[3];
            for (int j = 0; j < 3; j++)
            {
                double sum = 0.0;
                for (int i = 0; i < n; i++) sum += Q[i, j] * y[i];
                qty[j] = sum;
            }

            double[] coeffs = new double[3];
            for (int i = 2; i >= 0; i--)
            {
                double s = qty[i];
                for (int j = i + 1; j < 3; j++) s -= R[i, j] * coeffs[j];
                if (Math.Abs(R[i, i]) < 1e-15) return null;
                coeffs[i] = s / R[i, i];
            }

            return coeffs;
        }

        static List<string> FindExpressions(double[] xAxis, double[] yAxis, int maxDepth, int minConst, int maxConst, double tolerance)
        {
            var results = new List<string>();

            if (xAxis.Length == yAxis.Length && xAxis.Zip(yAxis, (x, y) => Math.Abs(x - y) <= tolerance).All(b => b))
            {
                results.Add("x");
                return results;
            }

            var ops = Enum.GetValues(typeof(ArithmeticActionType)).Cast<ArithmeticActionType>().ToArray();

            var constants = new List<double>();
            for (int c = minConst; c <= maxConst; c++) constants.Add(c);

            for (int depth = 1; depth <= maxDepth; depth++)
            {
                GenerateSequences(depth, ops, constants, seq =>
                {
                    if (IsMatchSequence(seq, xAxis, yAxis, tolerance, out string? repr))
                    {
                        results.Add(repr ?? "");
                    }
                });
            }

            return results.Distinct().ToList();
        }

        record Step(ArithmeticActionType Op, double Constant);

        static void GenerateSequences(int depth, ArithmeticActionType[] ops, List<double> constants, Action<List<Step>> onSequence)
        {
            var seq = new List<Step>(depth);

            void Recurse(int i)
            {
                if (i == depth)
                {
                    onSequence(new List<Step>(seq));
                    return;
                }

                foreach (var op in ops)
                {
                    foreach (var c in constants)
                    {
                        seq.Add(new Step(op, c));
                        Recurse(i + 1);
                        seq.RemoveAt(seq.Count - 1);
                    }
                }
            }

            Recurse(0);
        }

        static bool IsMatchSequence(List<Step> seq, double[] xAxis, double[] yAxis, double tol, out string? representation)
        {
            representation = null;

            var action = new ArithmeticAction();

            var results = new double[xAxis.Length];
            for (int i = 0; i < xAxis.Length; i++)
            {
                double val = xAxis[i];
                bool ok = true;
                foreach (var s in seq)
                {
                    if (s.Op == ArithmeticActionType.Div && Math.Abs(s.Constant) < double.Epsilon)
                    {
                        ok = false;
                        break;
                    }

                    action.Type = s.Op;
                    try
                    {
                        val = action.Perform(val, s.Constant);
                    }
                    catch
                    {
                        ok = false;
                        break;
                    }

                    if (double.IsNaN(val) || double.IsInfinity(val))
                    {
                        ok = false;
                        break;
                    }
                }

                if (!ok) return false;
                results[i] = val;
            }

            for (int i = 0; i < results.Length; i++)
            {
                if (Math.Abs(results[i] - yAxis[i]) > tol) return false;
            }

            string repr = "x";
            foreach (var s in seq)
            {
                string constStr = s.Constant % 1 == 0 ? ((int)s.Constant).ToString() : s.Constant.ToString();
                switch (s.Op)
                {
                    case ArithmeticActionType.Add: repr = $"({repr} + {constStr})"; break;
                    case ArithmeticActionType.Sub: repr = $"({repr} - {constStr})"; break;
                    case ArithmeticActionType.Mul: repr = $"({repr} * {constStr})"; break;
                    case ArithmeticActionType.Div: repr = $"({repr} / {constStr})"; break;
                    case ArithmeticActionType.Pow: repr = $"({repr} ^ {constStr})"; break;
                }
            }

            representation = repr;
            return true;
        }
    }
}
