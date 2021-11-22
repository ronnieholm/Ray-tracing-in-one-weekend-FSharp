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

## Benchmarking

The F# implementation is based on newer version of the book than the C# version.
So the benchmark isn't comparable:

.NET 5 (Linux)

```console
$ dotnet build --configuration release
$ dotnet run --configuration release RayTracingInOneWeekend.fsproj benchmark
Render time: 00:01:33.0886140
```

## References

- [3D Programming Fundamentals [3D Perspective Projection] Tutorial 4](https://www.youtube.com/watch?v=UgM6mIQfGDA&list=PLqCJpWy5Fohe8ucwhksiv9hTF5sfid8lA&index=5)