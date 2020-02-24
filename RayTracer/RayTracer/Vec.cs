using System;
using System.Collections.Generic;
using System.Text;

namespace RayTracer
{
    class Vec
    {
        public double X { get; set; }
        public double Y { get; set; }

        public double Z { get; set; }

        public Vec(double x=0, double y=0, double z=0)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public static Vec Normalize(Vec vec) { return vec.GetNormal(); }

        public static Vec Cross(Vec right, Vec left) { return right.CrossWith(left); }

        public static double Dot(Vec right, Vec left) { return right.DotWith(left); }

        public static Vec Multiply(Vec right, Vec left) { return right * left; }

        public static Vec Multiply(Vec right, double left) { return right * left; }

        public static Vec Add(Vec right, Vec left) { return right + left; }

        public static Vec Subtract(Vec right, Vec left) { return right - left; }
    
        public static Vec operator+(Vec right, Vec left) { return new Vec(right.X + left.X, right.Y + left.Y, right.Z + left.Z); }
        public static Vec operator-(Vec right, Vec left) { return new Vec(right.X - left.X, right.Y - left.Y, right.Z - left.Z); }

        public static Vec operator *(Vec right, Vec left) { return new Vec(right.X * left.X, right.Y * left.Y, right.Z * left.Z); }

        public static Vec operator *(Vec right, double left) { return new Vec(right.X * left, right.Y * left, right.Z * left); }

        public Vec GetNormal() { return this * (1 / Math.Sqrt(X * X + Y * Y + Z * Z)); }

        public double DotWith(Vec b) { return X * b.X + Y * b.Y + Z * b.Z; }

        public Vec CrossWith(Vec b) { return new Vec(Y * b.Z - Z * b.Y, Z * b.X - X * b.Z, X * b.Y - Y * b.X); }
    }
}
