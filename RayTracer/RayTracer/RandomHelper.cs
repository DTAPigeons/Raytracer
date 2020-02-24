using System;
using System.Collections.Generic;
using System.Text;

namespace RayTracer
{
    class RandomHelper
    {
        private Random random;
        private int seed;

        public RandomHelper(int seed)
        {
            random = new Random(seed);
            this.seed = seed;
        }

        public double NextDouble(int seed)
        {

            if(this.seed != seed)
            {
                random = new Random(seed);
                this.seed = seed;
            }

            return random.NextDouble();
        }
    }
}
