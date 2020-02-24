using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

namespace RayTracer
{
    class Scene
    {
        private int with = 1024;
        private int height = 768;
        private Ray camera = new Ray(new Vec(50.0, 52.0, 295.6), Vec.Normalize(new Vec(0f, -0.042612, -1)));

        private RandomHelper random = new RandomHelper(0);

        public Scene() { }

        public Scene(int with, int height, Ray camera, Sphere[] spheres)
        {
            this.with = with;
            this.height = height;
            this.camera = camera;
            this.spheres = spheres;
        }

        public Sphere[] spheres = {
            new Sphere(1e5, new Vec(1e5+1, 40.8, 81.6), new Vec(), new Vec(.75,.25,.25), Sphere.ReflectionType.DIFFUSE),             //Left
            new Sphere(1e5, new Vec(-1e5+99, 40.8, 81.6), new Vec(), new Vec(.25,.25,.75), Sphere.ReflectionType.DIFFUSE),           //Right
            new Sphere(1e5, new Vec(50, 40.8, 1e5),new Vec(), new Vec(.75,.75,.75), Sphere.ReflectionType.DIFFUSE),                 //Back
            new Sphere(1e5, new Vec(50, 40.8, -1e5+170),new Vec(),new Vec(), Sphere.ReflectionType.DIFFUSE),                          //Front
            new Sphere(1e5, new Vec(50, 1e5, 81.6),new Vec(), new Vec(.75,.75,.75), Sphere.ReflectionType.DIFFUSE),                 //Botom
            new Sphere(1e5, new Vec(50, -1e5+81.6, 81.6),new Vec(), new Vec(.75,.75,.75), Sphere.ReflectionType.DIFFUSE),          //Top
            new Sphere(16.5, new Vec(27, 16.5, 47),new Vec(), Vec.Multiply(new Vec(1,1,1), .999), Sphere.ReflectionType.SPECULAR),   //Mirror
            new Sphere(16.5, new Vec(73, 16.5, 78),new Vec(), Vec.Multiply(new Vec(1,1,1), .999), Sphere.ReflectionType.REFRACTIVE), //Glass
            new Sphere(600, new Vec(50, 681.6-.27, 81.6), new Vec(12,12,12),new Vec(), Sphere.ReflectionType.DIFFUSE),                //Lite
        };

        public Vec[] TraceScene(int samps)
        {
            int samples = samps / 4;
            Vec xDirectionIncrement = new Vec(((double)with * .5135 / (double)height),0,0);
            Vec yDirectionIncrement = Vec.Multiply(Vec.Normalize(Vec.Cross(xDirectionIncrement, camera.Direction)), .5135);
            Vec colorVector =new Vec();
            Vec[] image = new Vec[with * height];

            for(int i =0; i<image.Length; i++)
            {
                image[i] =new Vec();
            }

            for(int y=0; y<height; y++)
            {
                Console.WriteLine("Tracing: " + y);
                int seed = (y * y * y);
                for(int x =0; x<with; x++)
                {
                    if(y==13 && x >= 500)
                    {
                        int h = 0;
                    }
                    for(int sy=0, i= (height - y - 1) * with + x; sy<2; sy++)
                    {
                        for (int sx = 0; sx < 2; sx++, colorVector =new Vec())
                        {
                            for (int s = 0; s < samples; s++)
                            {
                                double r1 = 2.0 * random.NextDouble(seed);
                                double dx = r1 < 1.0 ? Math.Sqrt(r1) - 1.0 : 1.0 - Math.Sqrt(2.0 - r1);
                                double r2 = 2.0 * random.NextDouble(seed);
                                double dy = r2 < 1.0 ? Math.Sqrt(r2) - 1.0 : 1.0 - Math.Sqrt(2.0 - r2);

                                Vec part1 = Vec.Multiply(xDirectionIncrement, (((sx + .5 + dx) / 2.0 + x) / with - .5));
                                Vec part2 = Vec.Multiply(yDirectionIncrement, (((sy + .5 + dy) / 2.0 + y) / height - .5));

                                Vec radianceDirection = Vec.Add(part1, part2);
                                radianceDirection = Vec.Add(radianceDirection, camera.Direction);

                                Vec radianceOrigin = Vec.Multiply(radianceDirection, 140);
                                radianceOrigin = Vec.Add(camera.Origin, radianceOrigin);

                                Ray radianceRay = new Ray(radianceOrigin, Vec.Normalize(radianceDirection));

                                Vec radiance = ComputeRadiance(radianceRay, 0, seed);
                                radiance = Vec.Multiply(radiance, 1.0 / (double)samples);


                                colorVector = Vec.Add(colorVector, radiance);
                            }



                            Vec clampedColor = new Vec(ClampDoulbe(colorVector.X), ClampDoulbe(colorVector.Y), ClampDoulbe(colorVector.Z));
                            clampedColor = Vec.Multiply(clampedColor, .25);

                            image[i] = Vec.Add(image[i], clampedColor);
                            if(double.IsNaN(image[i].X) || double.IsNaN(image[i].Y) || double.IsNaN(image[i].Z))
                            {
                                throw new Exception();
                            }
                        }
                    }
                }
            }

            return image;
        }

        private Vec ComputeRadiance(Ray ray, int depth, int seed)
        {
            RayIntersectionData data = GetRayIntersectionData(ray);

            
            if (!data.Intersects || depth > 300) {
                return new Vec();
            }
            
            Sphere sphere = spheres[data.ClosestSpereId];

            double distanceToIntersection = data.Distance;
            int id = data.ClosestSpereId;

            Vec intersectionPoint = ray.Origin + Vec.Multiply(ray.Direction, distanceToIntersection);
          //  intersectionPoint = intersectionPoint + Vec.Multiply(Vec.Normalize(intersectionPoint), 0.05);
            Vec sphereNormal = Vec.Normalize(Vec.Subtract(intersectionPoint, sphere.Position));
            Vec surfaceNormal = Vec.Dot(sphereNormal, ray.Direction) < 0 ? sphereNormal : Vec.Multiply(sphereNormal, -1.0);
            Vec color = sphere.Color;


            double maxRefl = color.X > color.Y && color.X > color.Z ? color.X : color.Y > color.Z ? color.Y : color.Z;

            depth++;

            if (depth > 5 || maxRefl==0.0)
            {
                double randNum = random.NextDouble(seed);
                if (randNum < maxRefl ) {
                    color = Vec.Multiply(color, 1.0 / maxRefl);
                }
                else {
                    return sphere.Emission;
                }
            }

            if(sphere.Reflection == Sphere.ReflectionType.DIFFUSE)
            {
                double angleRand = random.NextDouble(seed) *2*Math.PI;
                double distanceRand = random.NextDouble(seed);
                double distanceRandSqtr = Math.Sqrt(distanceRand);

                Vec w = surfaceNormal;
                Vec u = Vec.Normalize(Vec.Cross(Math.Abs(w.X) > .1 ? new Vec(0, 1, 0) : new Vec(1, 0, 0), w));
                Vec v = Vec.Cross(w, u);

                Vec ref1 = Vec.Multiply(u, Math.Cos(angleRand));
                ref1 = Vec.Multiply(ref1, distanceRandSqtr);
                Vec ref2 = Vec.Multiply(v, Math.Sin(angleRand));
                ref2 = Vec.Multiply(ref2, distanceRandSqtr);
                Vec ref3 = Vec.Multiply(w, Math.Sqrt(1 - distanceRand));
                Vec ref4 = Vec.Add(ref1, ref2);
                ref4 = Vec.Add(ref4, ref3);

                Vec reflectionRayRand = Vec.Normalize(ref4);

                Vec nextRadiance = ComputeRadiance(new Ray(intersectionPoint, reflectionRayRand), depth, seed);

                Vec result = Vec.Multiply(color, nextRadiance);
                result = Vec.Add(sphere.Emission, result);

                if (double.IsNaN(result.X) || double.IsNaN(result.Y) || double.IsNaN(result.Z))
                {
                    throw new Exception();
                } 

                return result;
            }
            else if(sphere.Reflection == Sphere.ReflectionType.SPECULAR)
            {
                //FUCK UP POINT
                Vec ray1 = Vec.Multiply(sphereNormal, 2);
                double dot1 = Vec.Dot(sphereNormal, ray.Direction);
                ray1 = Vec.Multiply(ray1, dot1);
                Vec ray2 = Vec.Subtract(ray.Direction, ray1);

                Vec nextRadiance = ComputeRadiance(new Ray(intersectionPoint, ray2), depth, seed);

                Vec result = Vec.Multiply(color, nextRadiance);
                result = Vec.Add(sphere.Emission, result);

                if (double.IsNaN(result.X) || double.IsNaN(result.Y) || double.IsNaN(result.Z))
                {
                    throw new Exception();
                }

                return result;
            }
            else
            {
                //FUCK UP POINT
                Vec ray1 = Vec.Multiply(sphereNormal, 2f);
                double dot1 = Vec.Dot(sphereNormal, ray.Direction);
                ray1 = Vec.Multiply(ray1, dot1);
                Vec ray2 = Vec.Subtract(ray.Direction, ray1);

                Ray reflectionRay = new Ray(intersectionPoint, ray2);
                bool goesInto = Vec.Dot(sphereNormal, surfaceNormal)>0;
                double nc = 1.0;
                double nt = 1.5;
                double nnt = goesInto ? nc / nt : nt / nc;
                double ddn = Vec.Dot(ray.Direction, surfaceNormal);
                double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

                if (cos2t < 0)
                {
                    Vec nextRadiance = ComputeRadiance(reflectionRay, depth, seed);

                    Vec result = Vec.Multiply(color, nextRadiance);
                    result = Vec.Add(sphere.Emission, result);

                    if (double.IsNaN(result.X) || double.IsNaN(result.Y) || double.IsNaN(result.Z))
                    {
                        throw new Exception();
                    }

                    return result;
                }
                else
                {
                    Vec part1 = Vec.Multiply(ray.Direction, nnt);
                    double goesIntoMultiplier = goesInto ? 1.0 : -1.0;
                    Vec part2 = Vec.Multiply(sphereNormal, goesIntoMultiplier);
                    Vec part3 = Vec.Multiply(part2, (ddn * nnt + Math.Sqrt(cos2t)));
                    Vec part4 = Vec.Subtract(part1, part3);
                    Vec travelDirection = Vec.Normalize(part4);

                    double a = nt - nc;
                    double b = nt + nc;
                    double R0 = a * a / (b * b);
                    double c = 1.0 - (goesInto ? -ddn : Vec.Dot(travelDirection, sphereNormal));
                    double Re = R0 + (1 - R0) * c * c * c * c * c;
                    double Tr = 1.0 - Re;
                    double P = .25 + .5 * Re;
                    double RP = Re / P;
                    double TP = Tr / (1.0 - P);

                    Vec radianceUsed =new Vec();

                    if (depth > 2)
                    {
                        if(random.NextDouble(seed) < P) {
                            radianceUsed = ComputeRadiance(reflectionRay, depth, seed);
                            radianceUsed = Vec.Multiply(radianceUsed, RP);
                        }
                        else
                        {
                            radianceUsed = ComputeRadiance(new Ray(intersectionPoint, travelDirection), depth, seed);
                            radianceUsed = Vec.Multiply(radianceUsed, TP);
                        }
                    }
                    else
                    {
                        Vec nextRadiance1 = ComputeRadiance(reflectionRay, depth, seed);
                        nextRadiance1 = Vec.Multiply(nextRadiance1, Re);
                        Vec nextRadiance2 = ComputeRadiance(new Ray(intersectionPoint, travelDirection), depth, seed);
                        nextRadiance2 = Vec.Multiply(nextRadiance2, Tr);

                        radianceUsed = Vec.Add(nextRadiance1, nextRadiance2);
                    }

                    Vec result = Vec.Multiply(color, radianceUsed);
                    result = Vec.Add(sphere.Emission, result);

                    if (double.IsNaN(result.X) || double.IsNaN(result.Y) || double.IsNaN(result.Z))
                    {
                        throw new Exception();
                    }

                    return result;
                }

            }
        }

        public RayIntersectionData GetRayIntersectionData(Ray ray)
        {
            int numberOfSpheres = spheres.Length;
            int closestSpereId = -1;
            double inf = 1e20;
            double shortestDistance = inf;


            for(int i=numberOfSpheres-1; i>=0; i--)
            {
                double intersectionDistance = spheres[i].IntersectWithRay(ray);

                if(intersectionDistance!=0.0 && shortestDistance > intersectionDistance)
                {
                    shortestDistance = intersectionDistance;
                    closestSpereId = i;
                }
            }

            return new RayIntersectionData()
            {
                Ray = ray,
                Distance = shortestDistance,
                ClosestSpereId = closestSpereId,
                Intersects = shortestDistance < inf
            };
        }



        private double ClampDoulbe(double x)
        {
            return x < 0.0 ? 0.0 : x > 1.0 ? 1.0 : x;
        }
    }
}
