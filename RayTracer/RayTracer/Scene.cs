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
        private Ray camera = new Ray(new Vector3(50f, 52f, 295.6f), Vector3.Normalize(new Vector3(0f, -0.042612f, -1f)));

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
            new Sphere(1e5, new Vector3((float)1e5+1f, 40.8f, 81.6f), Vector3.Zero, new Vector3(.75f,.25f,.25f), Sphere.ReflectionType.DIFFUSE),             //Left
            new Sphere(1e5, new Vector3((float)-1e5+99f, 40.8f, 81.6f), Vector3.Zero, new Vector3(.25f,.25f,.75f), Sphere.ReflectionType.DIFFUSE),           //Right
            new Sphere(1e5, new Vector3(50f, 40.8f, (float)1e5), Vector3.Zero, new Vector3(.75f,.75f,.75f), Sphere.ReflectionType.DIFFUSE),                 //Back
            new Sphere(1e5, new Vector3(50f, 40.8f, (float)-1e5+170f), Vector3.Zero, Vector3.Zero, Sphere.ReflectionType.DIFFUSE),                          //Front
            new Sphere(1e5, new Vector3(50f, (float)1e5, 81.6f), Vector3.Zero, new Vector3(.75f,.75f,.75f), Sphere.ReflectionType.DIFFUSE),                 //Botom
            new Sphere(1e5, new Vector3(50f, (float)-1e5+81.6f, 81.6f), Vector3.Zero, new Vector3(.75f,.75f,.75f), Sphere.ReflectionType.DIFFUSE),          //Top
            new Sphere(16.5, new Vector3(27f, 16.5f, 47f), Vector3.Zero, Vector3.Multiply(new Vector3(1f,1f,1f), .999f), Sphere.ReflectionType.SPECULAR),   //Mirror
            new Sphere(16.5, new Vector3(73f, 16.5f, 78f), Vector3.Zero, Vector3.Multiply(new Vector3(1f,1f,1f), .999f), Sphere.ReflectionType.REFRACTIVE), //Glass
            new Sphere(600, new Vector3(50f, 681.6f-.27f, 81.6f), new Vector3(12f,12f,12f), Vector3.Zero, Sphere.ReflectionType.DIFFUSE),                //Lite
        };

        public Vector3[] TraceScene(int samps)
        {
            int samples = samps / 4;
            Vector3 xDirectionIncrement = new Vector3((float)((double)with * .5135 / (double)height),0f,0f);
            Vector3 yDirectionIncrement = Vector3.Multiply(Vector3.Normalize(Vector3.Cross(xDirectionIncrement, camera.Direction)), .5135f);
            Vector3 colorVector = Vector3.Zero;
            Vector3[] image = new Vector3[with * height];

            for(int i =0; i<image.Length; i++)
            {
                image[i] = Vector3.Zero;
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
                        for (int sx = 0; sx < 2; sx++, colorVector = Vector3.Zero)
                        {
                            for (int s = 0; s < samples; s++)
                            {
                                double r1 = 2.0 * random.NextDouble(seed);
                                double dx = r1 < 1.0 ? Math.Sqrt(r1) - 1.0 : 1.0 - Math.Sqrt(2.0 - r1);
                                double r2 = 2.0 * random.NextDouble(seed);
                                double dy = r2 < 1.0 ? Math.Sqrt(r2) - 1.0 : 1.0 - Math.Sqrt(2.0 - r2);

                                Vector3 part1 = Vector3.Multiply(xDirectionIncrement, (float)(((sx + .5 + dx) / 2.0 + x) / with - .5));
                                Vector3 part2 = Vector3.Multiply(yDirectionIncrement, (float)(((sy + .5 + dy) / 2.0 + y) / height - .5));

                                Vector3 radianceDirection = Vector3.Add(part1, part2);
                                radianceDirection = Vector3.Add(radianceDirection, camera.Direction);

                                Vector3 radianceOrigin = Vector3.Multiply(radianceDirection, 140f);
                                radianceOrigin = Vector3.Add(camera.Origin, radianceOrigin);

                                Ray radianceRay = new Ray(radianceOrigin, Vector3.Normalize(radianceDirection));

                                Vector3 radiance = ComputeRadiance(radianceRay, 0, seed);
                                radiance = Vector3.Multiply(radiance, 1f / samples);


                                colorVector = Vector3.Add(colorVector, radiance);
                            }



                            Vector3 clampedColor = new Vector3((float)ClampDoulbe(colorVector.X), (float)ClampDoulbe(colorVector.Y), (float)ClampDoulbe(colorVector.Z));
                            clampedColor = Vector3.Multiply(clampedColor, .25f);

                            image[i] = Vector3.Add(image[i], clampedColor);
                            if(float.IsNaN(image[i].X) || float.IsNaN(image[i].Y) || float.IsNaN(image[i].Z))
                            {
                                throw new Exception();
                            }
                        }
                    }
                }
            }

            return image;
        }

        private Vector3 ComputeRadiance(Ray ray, int depth, int seed)
        {
            RayIntersectionData data = GetRayIntersectionData(ray);

            
            if (!data.Intersects || depth > 300) {
                return Vector3.Zero;
            }
            
            Sphere sphere = spheres[data.ClosestSpereId];

            double distanceToIntersection = data.Distance;
            int id = data.ClosestSpereId;

            Vector3 intersectionPoint = ray.Origin + Vector3.Multiply(ray.Direction, (float)distanceToIntersection);
          //  intersectionPoint = intersectionPoint + Vector3.Multiply(Vector3.Normalize(intersectionPoint), 0.05f);
            Vector3 sphereNormal = Vector3.Normalize(Vector3.Subtract(intersectionPoint, sphere.Position));
            Vector3 surfaceNormal = Vector3.Dot(sphereNormal, ray.Direction) < 0f ? sphereNormal : Vector3.Multiply(sphereNormal, -1f);
            Vector3 color = sphere.Color;


            double maxRefl = color.X > color.Y && color.X > color.Z ? color.X : color.Y > color.Z ? color.Y : color.Z;

            depth++;

            if (depth > 5 || maxRefl==0.0)
            {
                double randNum = random.NextDouble(seed);
                if (randNum < maxRefl ) {
                    color = Vector3.Multiply(color, 1f / (float)maxRefl);
                }
                else {
                    return sphere.Emission;
                }
            }

            if(sphere.Reflection == Sphere.ReflectionType.DIFFUSE)
            {
                double angleRand = random.NextDouble(seed) *2f*Math.PI;
                double distanceRand = random.NextDouble(seed);
                double distanceRandSqtr = Math.Sqrt(distanceRand);

                Vector3 w = surfaceNormal;
                Vector3 u = Vector3.Normalize(Vector3.Cross(Math.Abs(w.X) > .1 ? new Vector3(0f, 1f, 0f) : new Vector3(1f, 0f, 0f), w));
                Vector3 v = Vector3.Cross(w, u);

                Vector3 ref1 = Vector3.Multiply(u, (float)Math.Cos(angleRand));
                ref1 = Vector3.Multiply(ref1, (float)distanceRandSqtr);
                Vector3 ref2 = Vector3.Multiply(v, (float)Math.Sin(angleRand));
                ref2 = Vector3.Multiply(ref2, (float)distanceRandSqtr);
                Vector3 ref3 = Vector3.Multiply(w, (float)Math.Sqrt(1 - distanceRand));
                Vector3 ref4 = Vector3.Add(ref1, ref2);
                ref4 = Vector3.Add(ref4, ref3);

                Vector3 reflectionRayRand = Vector3.Normalize(ref4);

                Vector3 nextRadiance = ComputeRadiance(new Ray(intersectionPoint, reflectionRayRand), depth, seed);

                Vector3 result = Vector3.Multiply(color, nextRadiance);
                result = Vector3.Add(sphere.Emission, result);

                if (float.IsNaN(result.X) || float.IsNaN(result.Y) || float.IsNaN(result.Z))
                {
                    throw new Exception();
                } 

                return result;
            }
            else if(sphere.Reflection == Sphere.ReflectionType.SPECULAR)
            {
                //FUCK UP POINT
                Vector3 ray1 = Vector3.Multiply(sphereNormal, 2f);
                float dot1 = Vector3.Dot(sphereNormal, ray.Direction);
                ray1 = Vector3.Multiply(ray1, dot1);
                Vector3 ray2 = Vector3.Subtract(ray.Direction, ray1);

                Vector3 nextRadiance = ComputeRadiance(new Ray(intersectionPoint, ray2), depth, seed);

                Vector3 result = Vector3.Multiply(color, nextRadiance);
                result = Vector3.Add(sphere.Emission, result);

                if (float.IsNaN(result.X) || float.IsNaN(result.Y) || float.IsNaN(result.Z))
                {
                    throw new Exception();
                }

                return result;
            }
            else
            {
                //FUCK UP POINT
                Vector3 ray1 = Vector3.Multiply(sphereNormal, 2f);
                float dot1 = Vector3.Dot(sphereNormal, ray.Direction);
                ray1 = Vector3.Multiply(ray1, dot1);
                Vector3 ray2 = Vector3.Subtract(ray.Direction, ray1);

                Ray reflectionRay = new Ray(intersectionPoint, ray2);
                bool goesInto = Vector3.Dot(sphereNormal, surfaceNormal)>0;
                double nc = 1.0;
                double nt = 1.5;
                double nnt = goesInto ? nc / nt : nt / nc;
                double ddn = Vector3.Dot(ray.Direction, surfaceNormal);
                double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

                if (cos2t < 0)
                {
                    Vector3 nextRadiance = ComputeRadiance(reflectionRay, depth, seed);

                    Vector3 result = Vector3.Multiply(color, nextRadiance);
                    result = Vector3.Add(sphere.Emission, result);

                    if (float.IsNaN(result.X) || float.IsNaN(result.Y) || float.IsNaN(result.Z))
                    {
                        throw new Exception();
                    }

                    return result;
                }
                else
                {
                    Vector3 part1 = Vector3.Multiply(ray.Direction, (float)nnt);
                    float goesIntoMultiplier = goesInto ? 1.0f : -1.0f;
                    Vector3 part2 = Vector3.Multiply(sphereNormal, goesIntoMultiplier);
                    Vector3 part3 = Vector3.Multiply(part2, (float)(ddn * nnt + Math.Sqrt(cos2t)));
                    Vector3 part4 = Vector3.Subtract(part1, part3);
                    Vector3 travelDirection = Vector3.Normalize(part4);

                    double a = nt - nc;
                    double b = nt + nc;
                    double R0 = a * a / (b * b);
                    double c = 1.0 - (goesInto ? -ddn : Vector3.Dot(travelDirection, sphereNormal));
                    double Re = R0 + (1 - R0) * c * c * c * c * c;
                    double Tr = 1.0 - Re;
                    double P = .25 + .5 * Re;
                    double RP = Re / P;
                    double TP = Tr / (1.0 - P);

                    Vector3 radianceUsed = Vector3.Zero;

                    if (depth > 2)
                    {
                        if(random.NextDouble(seed) < P) {
                            radianceUsed = ComputeRadiance(reflectionRay, depth, seed);
                            radianceUsed = Vector3.Multiply(radianceUsed, (float)RP);
                        }
                        else
                        {
                            radianceUsed = ComputeRadiance(new Ray(intersectionPoint, travelDirection), depth, seed);
                            radianceUsed = Vector3.Multiply(radianceUsed, (float)TP);
                        }
                    }
                    else
                    {
                        Vector3 nextRadiance1 = ComputeRadiance(reflectionRay, depth, seed);
                        nextRadiance1 = Vector3.Multiply(nextRadiance1, (float)Re);
                        Vector3 nextRadiance2 = ComputeRadiance(new Ray(intersectionPoint, travelDirection), depth, seed);
                        nextRadiance2 = Vector3.Multiply(nextRadiance2, (float)Tr);

                        radianceUsed = Vector3.Add(nextRadiance1, nextRadiance2);
                    }

                    Vector3 result = Vector3.Multiply(color, radianceUsed);
                    result = Vector3.Add(sphere.Emission, result);

                    if (float.IsNaN(result.X) || float.IsNaN(result.Y) || float.IsNaN(result.Z))
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
