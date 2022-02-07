# TODO

- Make rayColor function non-recursive/iterative. Current implementation isn't tail 
  recursive and performance suffers.
- Use https://sharplab.io to compare source code
- Use Benchmark.NET (and memory benchmarks) on per pixel operations in F#/C#
- Try out [<MethodImpl(MethodImplOptions.AggressiveOptimization>]