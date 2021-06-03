# Ray-tracing-in-one-weekend-FSharp

Based on [Ray Tracing in One Weekend](https://raytracing.github.io/books/RayTracingInOneWeekend.html), Version 3.2.3,
2020-12-07.

See also [C# implementation](https://github.com/ronnieholm/Ray-tracing-in-one-weekend-CSharp). 

![Random scene](Random-scene.png)

## Getting started

    $ git clone https://github.com/ronnieholm/Ray-tracing-in-one-weekend-FSharp.git
    $ dotnet build
    $ dotnet run --configuration release RayTracingInOneWeekend.fsproj > out.ppm
    $ display out.ppm