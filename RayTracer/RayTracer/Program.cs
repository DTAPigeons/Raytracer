using System;
using System.IO;
using System.Numerics;
using System.Text;

namespace RayTracer
{
    class Program
    {
        private static int GetDisplayInt(double x)
        {
            return (int)(Math.Pow(ClampDoulbe(x), 1.0 / 2.2) * 255.0 + .5);
        }

        private static double ClampDoulbe(double x)
        {
            return x < 0.0 ? 0.0 : x > 1.0 ? 1.0 : x;
        }

        static void Main(string[] args)
        {
            string path = @"result.ppm";

            Scene scene = new Scene();

            Vector3[] image = scene.TraceScene(4);

            StringBuilder sb = new StringBuilder();

            foreach(Vector3 pixel in image)
            {
                int x = GetDisplayInt(pixel.X);
                int y = GetDisplayInt(pixel.Y);
                int z = GetDisplayInt(pixel.Z);
                sb.Append(x + " " + y + " " + z + " ");

            }

            using (StreamWriter sw = File.CreateText(path))
            {
                sw.WriteLine("P3");
                sw.WriteLine("1024 768");
                sw.WriteLine("255");
                sw.WriteLine(sb.ToString());
            }
            
        }
    }
}
