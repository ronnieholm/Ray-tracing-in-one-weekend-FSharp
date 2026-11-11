# Ray-tracing-in-one-weekend-FSharp

Based on [Ray Tracing in One Weekend](https://raytracing.github.io/books/RayTracingInOneWeekend.html), Version 3.2.3, 2020-12-07.

See also [C# implementation](https://github.com/ronnieholm/Ray-tracing-in-one-weekend-CSharp).

![Random scene](Random-scene.png)

## Getting started

    $ git clone https://github.com/ronnieholm/Ray-tracing-in-one-weekend-FSharp.git
    $ dotnet run --configuration release RayTracingInOneWeekend.fsproj > out.ppm
    $ display out.ppm

## Benchmarking

The F# implementation is based on newer version of the book than the C# version.
So benchmarks aren't comparable:

Benchmarks run on AMD Ryzen 9 9950X 16-Core processor.

.NET 11 (Linux)

    $ /bin/time dotnet run --configuration release RayTracingInOneWeekend.fsproj benchmark
    Render time: 00:00:30.3248245
    30.85user 0.10system 0:30.77elapsed 100%CPU (0avgtext+0avgdata 187224maxresident)k
    0inputs+32outputs (0major+48837minor)pagefaults 0swaps

.NET 10 (Linux)

    $ /bin/time dotnet run --configuration release RayTracingInOneWeekend.fsproj benchmark
    Render time: 00:00:18.9031988
    19.88user 0.13system 0:19.59elapsed 102%CPU (0avgtext+0avgdata 196824maxresident)k
    0inputs+96outputs (0major+57150minor)pagefaults 0swaps

.NET 9 (Linux)

    $ /bin/time dotnet run --configuration release RayTracingInOneWeekend.fsproj benchmark
    Render time: 00:00:16.9815730
    17.73user 0.13system 0:17.60elapsed 101%CPU (0avgtext+0avgdata 173852maxresident)k
    0inputs+56outputs (0major+51335minor)pagefaults 0swaps

## References

- [C++ reference implementation](https://github.com/RayTracing/raytracing.github.io)
- [3D Programming Fundamentals [3D Perspective Projection] Tutorial 4](https://www.youtube.com/watch?v=UgM6mIQfGDA&list=PLqCJpWy5Fohe8ucwhksiv9hTF5sfid8lA&index=5)