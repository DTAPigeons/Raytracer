using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

namespace RayTracer
{

    class Sphere
    {
        private double radius;
        private Vec position;
        private Vec emission;
        private Vec color;
        private ReflectionType reflection;

        public Sphere(double radius, Vec position, Vec emission, Vec color, ReflectionType reflection)
        {
            this.Radius = radius;
            this.Position = position;
            this.Emission = emission;
            this.Color = color;
            this.Reflection = reflection;
        }

        public double Radius { get => radius; set => radius = value; }
        public Vec Position { get => position; set => position = value; }
        public Vec Emission { get => emission; set => emission = value; }
        public Vec Color { get => color; set => color = value; }
        internal ReflectionType Reflection { get => reflection; set => reflection = value; }

        public double IntersectWithRay(Ray ray)
        {
            // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
            Vec op = Position - ray.Origin;
            double t;
            double epsilon = 1e-4;
            double b = Vec.Dot(op, ray.Direction); // 1/2 b from the quadratic equesion
            double detetminant = b*b-op.DotWith(op)+Radius*Radius;
            if (detetminant < 0)
            {
                return 0.0;
            }
            else
            {
                detetminant = Math.Sqrt(detetminant);
                t = (b - detetminant > epsilon) ? b - detetminant : b + detetminant;
                t = t > epsilon ? t : 0.0;
                return t;
            }
        }



        public enum ReflectionType {
            DIFFUSE,
            SPECULAR,
            REFRACTIVE 

        }
    }
}
