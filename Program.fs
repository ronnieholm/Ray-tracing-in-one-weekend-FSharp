open System
open System.Diagnostics
open System.Runtime.CompilerServices

[<AutoOpen>]
module Utils =
    // Module wraps .NET RNG to abstract random number generation and to easily switch it out for a faster generator.
    let rng = Random()

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    let randomDouble (): double =
        rng.NextDouble()

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    let inline randomBetween (min: double) (max: double): double =
        min + (max - min) * randomDouble()

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    let degreesToRadians (degrees: double): double =
        degrees * Math.PI / 180.

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    let clamp (x: double) (min: double) (max: double): double =
        if x < min then min elif x > max then max else x

[<IsReadOnly; Struct>]
type Vec3 =
    { X: double
      Y: double
      Z: double }
with
    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    member v.Length() =
        Math.Sqrt(v.X * v.X + v.Y * v.Y + v.Z * v.Z)

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    member v.LengthSquared() =
        v.X * v.X + v.Y * v.Y + v.Z * v.Z

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    override v.ToString() =
        $"{v.X} {v.Y} {v.Z}"   

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member create x y z =
        { X = x; Y = y; Z = z }
        
    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member (~-) (v: Vec3) =
        { X = -v.X; Y = -v.Y; Z = -v.Z }

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member (+) (v: Vec3, u: Vec3) =
        { X = v.X + u.X; Y = v.Y + u.Y; Z = v.Z + u.Z }

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member (-) (v: Vec3, u: Vec3) =
        { X = v.X - u.X; Y = v.Y - u.Y; Z = v.Z - u.Z }

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member (*) (v: Vec3, u: Vec3) =
        { X = v.X * u.X; Y = v.Y * u.Y; Z = v.Z * u.Z }

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member (*) (v: Vec3, t: double) =
        { X = v.X * t; Y = v.Y * t; Z = v.Z * t }

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member (*) (t: double, v: Vec3) =
        { X = v.X * t; Y = v.Y * t; Z = v.Z * t }
        
    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member (/) (v: Vec3, t: double) =
        { X = v.X / t; Y = v.Y / t; Z = v.Z / t }

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member UnitVector (v: Vec3): Vec3 =
        v / v.Length()

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member Dot (v: Vec3) (u: Vec3) =
        v.X * u.X + v.Y * u.Y + v.Z * u.Z

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member Cross (v: Vec3) (u: Vec3) =
        { X = v.Y * u.Z - v.Z * u.Y
          Y = v.Z * u.X - v.X * u.Z
          Z = v.X * u.Y - v.Y * u.X }

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member Random(): Vec3 =
        { X = randomDouble(); Y = randomDouble(); Z = randomDouble() }

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member RandomBetween (min: double) (max: double): Vec3 =
        { X = randomBetween min max
          Y = randomBetween min max
          Z = randomBetween min max }

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member RandomInUnitSphere(): Vec3 =
        seq {
            while true do
                let p = Vec3.RandomBetween -1. 1.
                if p.LengthSquared() < 1.
                then yield p
        } |> Seq.take 1 |> Seq.head

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member RandomInUnitDisk() =        
        seq {
            while true do
                let p = Vec3.create (randomBetween -1. 1.) (randomBetween -1. 1.) 0.
                if p.LengthSquared() < 1.
                then yield p               
        } |> Seq.take 1 |> Seq.head

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member RandomUnitVector() =
        Vec3.UnitVector (Vec3.RandomInUnitSphere())

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    member v.NearZero(): bool =
        // Returns true if the vector is close to zero in all dimensions.
        let s = 1e-8
        Math.Abs v.X < s && Math.Abs v.Y < s && Math.Abs v.Z < s

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member Reflect (v: Vec3) (unitVector: Vec3): Vec3 =
        v - 2. * (Vec3.Dot v unitVector) * unitVector

    [<MethodImpl(MethodImplOptions.AggressiveInlining)>]
    static member Refract (unitVector: Vec3) (normal: Vec3) (etaIOverEtaT: double) =
        let cosTheta = Math.Min(Vec3.Dot -unitVector normal, 1.)
        let rOutPerpendicular = etaIOverEtaT * (unitVector + cosTheta * normal)
        let rOutParallel = -Math.Sqrt(Math.Abs(1. - rOutPerpendicular.LengthSquared())) * normal
        rOutPerpendicular + rOutParallel       

type Point3 = Vec3 // 3D point
type Color = Vec3  // RGB color

let black = Color.create 0. 0. 0.
let white = Color.create 1. 1. 1.
let blue = Color.create 0.5 0.7 1.0

[<IsReadOnly; Struct>]
type Ray =
    { Origin: Point3
      Direction: Vec3 }
with 
    member r.At(t: double): Point3 =
        r.Origin + t * r.Direction

let writeColor (c: Color) (samplesPerPixel: int) =
    // Divide the color by the number of samples and gamma-correct for gamma = 2
    let scale = 1. / double(samplesPerPixel)
    let r = sqrt (scale * c.X)
    let g = sqrt (scale * c.Y)
    let b = sqrt (scale * c.Z)
    printfn $"{int(256. * clamp r 0. 0.999)} {int(256. * clamp g 0. 0.999)} {int(256. * clamp b 0. 0.999)}"    

let hitSphere (center: Point3) (radius: double) (r: Ray) =
    let oc = r.Origin - center
    let a = r.Direction.LengthSquared()
    let halfB = Vec3.Dot oc r.Direction
    let c = oc.LengthSquared() - (radius * radius)
    let discriminant = halfB * halfB - a * c
    if discriminant < 0.
    then -1.
    else (-halfB - sqrt discriminant) / a

[<IsReadOnly; Struct>]
type HitRecord =
    { P: Point3
      Normal: Vec3
      Material: Material
      T: double
      FrontFacing: bool }

and Material =
   | Lambertian of albedo: Color
   | Metal of albedo: Color * fuzziness: double
   | Dielectric of refractionIndex: double
with
    // Returns (* attenuation *) Color * (* scattered *) Ray
    member m.scatter (r: Ray) (hit: HitRecord): Option<Color * Ray> =
        match m with
        | Lambertian albedo ->
            let scatterDirection = hit.Normal + Vec3.RandomUnitVector()
            let scatterDirection' = 
                if scatterDirection.NearZero() 
                then hit.Normal else scatterDirection
            Some (albedo, { Origin = hit.P; Direction = scatterDirection' })
        | Metal(albedo, fuzziness) ->
            let reflected = Vec3.Reflect (Vec3.UnitVector r.Direction) hit.Normal
            let scattered = { Origin = hit.P; Direction = reflected + fuzziness * Vec3.RandomInUnitSphere() }            
            if Vec3.Dot scattered.Direction hit.Normal > 0.            
            then Some (albedo, scattered)
            else None
        | Dielectric refractionIndex ->
            let reflectance (cosine: double) =
                // Use Schlick's approximation for reflectance.
                let r0 = (1. - refractionIndex) / (1. + refractionIndex)
                let r0' = r0 * r0
                r0' + (1. + r0') * Math.Pow(1. - cosine, 5.)
            
            let attenuation' = white
            let refractionRatio = if hit.FrontFacing then double(1. / refractionIndex) else refractionIndex
            let unitDirection = Vec3.UnitVector r.Direction
            let cosTheta = Math.Min(Vec3.Dot -unitDirection hit.Normal, 1.)
            let sinTheta = Math.Sqrt(1. - cosTheta * cosTheta)
            let cannotRefract = refractionRatio * sinTheta > 1.
            let direction =
                if cannotRefract || reflectance cosTheta > randomDouble()
                then Vec3.Reflect unitDirection hit.Normal
                else Vec3.Refract unitDirection hit.Normal refractionRatio                         
            Some (attenuation', { Origin = hit.P; Direction = direction })

type Hittable =
    | Sphere of center: Point3 * radius: double * material: Material
with
    member h.Hit (r: Ray) (tMin: double) (tMax: double): HitRecord option =
        let setFaceNormal (outwardNormal: Vec3) =
            let frontFacing = Vec3.Dot r.Direction outwardNormal < 0.
            let normal = if frontFacing then outwardNormal else -outwardNormal
            (frontFacing, normal)

        match h with
        | Sphere(center, radius, material) ->
            let oc = r.Origin - center
            let a = r.Direction.LengthSquared()
            let halfB = Vec3.Dot oc r.Direction
            let c = oc.LengthSquared() - (radius * radius)
            let discriminant = halfB * halfB - a * c
            if discriminant < 0. then None
            else
                let sqrtD = sqrt discriminant
                let hit (root: double) =                   
                    let p =  r.At root 
                    let frontFacing, normal = setFaceNormal ((p - center) / radius)
                    Some { T = root; P = p; Normal = normal; FrontFacing = frontFacing; Material = material }
                
                // Find the nearest root that lies in the acceptable range
                let mutable root = (-halfB - sqrtD) / a 
                if root < tMin || tMax < root then
                    root <- (-halfB + sqrtD) / a
                    if root < tMin || tMax < root
                    then None
                    else hit root
                else hit root

let randomScene =
    [| let groundMaterial = Lambertian(Color.create 0.5 0.5 0.5)
       yield Sphere((Point3.create 0. -1000. 0.), 1000., groundMaterial)
       
       for a = -11 to 10 do
           for b = -11 to 10 do
               let chooseMaterial = randomDouble()
               let center = Point3.create (double(a) + 0.9 * randomDouble()) 0.2 (double(b) + 0.9 * randomDouble())
               
               if (center - Point3.create 4. 0.2 0.).Length() > 0.9 then
                   if chooseMaterial < 0.8 then
                       // Diffuse
                       let albedo = Color.Random() * Color.Random()
                       let material = Lambertian(albedo)
                       yield Sphere(center, 0.2, material)
                   elif chooseMaterial < 0.95 then
                       // Metal
                       let albedo = Color.RandomBetween 0.5 1.
                       let fuzziness = randomBetween 0. 0.5
                       let material = Metal(albedo, fuzziness)
                       yield Sphere(center, 0.2, material)
                   else
                       // Glass
                       let material = Dielectric(1.5)
                       yield Sphere(center, 0.2, material)
                       
       let material1 = Dielectric(1.5)
       yield Sphere((Point3.create 0. 1. 0.), 1., material1)
       let material2 = Lambertian(Color.create 0.4 0.2 0.1)
       yield Sphere((Point3.create -4. 1. 0.), 1., material2)
       let material3 = Metal((Color.create 0.7 0.6 0.5), 0.)
       yield Sphere((Point3.create 4. 1. 0.), 1., material3) |]

let hit (r: Ray) (tMin: double) (tMax: double) (hittables: Hittable array): HitRecord option =
    let mutable hitRecord = None
    let mutable hitAnything = false
    let mutable closestSoFar = tMax

    for object in hittables do
        match object.Hit r tMin closestSoFar with
        | Some h ->
            hitAnything <- true
            closestSoFar <- h.T
            hitRecord <- Some h
        | None -> ()

    hitRecord

let rec rayColor (r: Ray) (world: Hittable array) (depth: int): Color =
    // If we've exceeded the ray bounce limit, no more light is gathered
    if depth <= 0 then black
    else
        match hit r 0.001 (Double.PositiveInfinity) world with
        | Some h ->
            let scatter = h.Material.scatter r h
            match scatter with
            | Some (attenuation, scatter) ->
                attenuation * rayColor scatter world (depth - 1)
            | None -> black
        | None ->
            let unitDirection = Vec3.UnitVector(r.Direction)
            let t = 0.5 * (unitDirection.Y + 1.0)

            // Note: a + t(b - a) = (1 t)a + tb
            (1.0 - t) * white + t * blue    

[<IsReadOnly; Struct>]
type Camera =
    { Origin: Point3
      Horizontal: Vec3
      Vertical: Vec3
      LowerLeftCorner: Point3
      U: Vec3
      V: Vec3
      W: Vec3
      LensRadius: double }
with 
    static member create (lookFrom: Point3) (lookAt: Point3) (viewUp: Vec3) (verticalFieldOfView (* in degrees *): double) (aspectRatio: double) (aperture: double) (focusDistance: double): Camera =
        let theta = degreesToRadians verticalFieldOfView
        let h = Math.Tan(theta / 2.)
        let viewportHeight = 2. * h
        let viewportWidth = aspectRatio * viewportHeight      

        let w = Vec3.UnitVector (lookFrom - lookAt)
        let u = Vec3.UnitVector (Vec3.Cross viewUp w)
        let v = Vec3.Cross w u

        let origin = lookFrom
        let horizontal =  focusDistance * viewportWidth * u
        let vertical = focusDistance * viewportHeight * v
        let lowerLeftCorner = origin - horizontal / 2. - vertical / 2. - focusDistance * w
        let lensRadius = aperture / 2.
        { Origin = origin
          Horizontal = horizontal
          Vertical = vertical
          LowerLeftCorner = lowerLeftCorner
          U = u
          V = v
          W = w
          LensRadius = lensRadius }
    member c.GetRay (s: double) (t: double): Ray =
        let rd = c.LensRadius * Vec3.RandomInUnitDisk()
        let offset = c.U * rd.X + c.V * rd.Y
        { Origin = c.Origin + offset
          Direction = c.LowerLeftCorner + s * c.Horizontal + t * c.Vertical - c.Origin - offset }                
    
// Image
let aspectRatio = 3. / 2.
let imageWidth = 1200
let imageHeight = int(double(imageWidth) / aspectRatio)
let samplePerPixel = 500
let maxDepth = 50

// Camera
let lookFrom = (Point3.create 13. 2. 3.)
let lookAt = (Point3.create 0. 0. 0.)
let viewUp = (Point3.create 0. 1. 0.)
let distanceToFocus = 10.
let aperture = 0.1
let camera = Camera.create lookFrom lookAt viewUp 20. aspectRatio aperture distanceToFocus              

[<EntryPoint>]
let main _ =
    let sw = Stopwatch()
    sw.Start()
    printf $"P3\n{imageWidth} {imageHeight}\n255\n"
    
    for j = imageHeight - 1 downto 0 do
        eprintf $"\rScan lines remaining: {j} "
        for i = 0 to imageWidth - 1 do
            let mutable pixelColor = black
            for _ = 0 to samplePerPixel - 1 do
                let u = (double(i) + randomDouble()) / double(imageWidth - 1)
                let v = (double(j) + randomDouble()) / double(imageHeight - 1)
                let r = camera.GetRay u v
                pixelColor <- pixelColor + rayColor r randomScene maxDepth
            writeColor pixelColor samplePerPixel

    eprintfn "\nDone."
    eprintfn $"Render time: {sw.Elapsed}"
    0