using System;
using System.Collections.Generic;
using System.Text;

namespace RayTracer
{
    class RayIntersectionData
    {
        public Ray Ray { get; set; }
        public double Distance { get; set; }

        public int ClosestSpereId { get; set; }

        public bool Intersects { get; set; }
    }
}
