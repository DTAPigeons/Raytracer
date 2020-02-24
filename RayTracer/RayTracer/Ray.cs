using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

namespace RayTracer
{
    class Ray
    {
        public Vec Origin { get; set; }
        public Vec Direction { get; set; }

        public Ray(Vec origin, Vec direction)
        {
            Origin = origin;
            Direction = direction;
        }
    }
}
