#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using SparkCL;
using MathShards.Matrices;
using Silk.NET.OpenCL;
using MathShards.Matrices.Types;

namespace Main.Tests.HalvesTests;

static class Common
{
    public static void HalfMultiplies(IHalves matrix, Real[] b) {
        var vec = new Real[b.Length];
        var vec2 = new Real[b.Length];
        // L
        b.AsSpan().CopyTo(vec);

        matrix.InvLMul(vec);
        matrix.LMul(vec, vec2);

        var err = vec2
            .Zip(b)
            .Sum(a =>
            {
                return Real.Abs(a.First - a.Second);
            });

        Assert.That(err, Is.EqualTo(0).Within(1e-12));

        // U
        b.AsSpan().CopyTo(vec);

        matrix.InvUMul(vec);
        matrix.UMul(vec, vec2);

        err = vec2
            .Zip(b)
            .Sum(a =>
            {
                return Real.Abs(a.First - a.Second);
            });

        Assert.That(err, Is.EqualTo(0).Within(1e-12));
    }
    
    public static void HalfMultipliesOpenCL(IHalves matrix, Real[] b) {

        var compMatrix = matrix.GetComputeMatrix() as SparkAlgos.Types.IHalves;
        var compB = new ComputeBuffer<Real>(b, BufferFlags.OnHostAndDevice);
        
        var comp1 = new ComputeBuffer<Real>(b.Length, BufferFlags.OnHostAndDevice);
        var vec1 = new Real[b.Length];
        
        // LMul
        compMatrix!.LMul(compB, comp1);
        matrix.LMul(b, vec1);
        comp1.ToHost();
        {
            // TODO: в прошлый раз я забыл сделать ToHost. Это знак,
            // что API непонятный
            using var accessor = comp1.MapHost(MapFlags.Read);
            for (int i = 0; i < b.Length; i++)
            {
                Assert.That(accessor[i], Is.EqualTo(vec1[i]).Within(1e-12));
            } 
        }

        // InvLMul
        compB.CopyDeviceTo(comp1);
        b.CopyTo(vec1);

        compMatrix.InvLMul(comp1);
        matrix.InvLMul(vec1);
        comp1.ToHost();
        {
            using var accessor = comp1.MapHost(MapFlags.Read);
            for (int i = 0; i < b.Length; i++)
            {
                Assert.That(accessor[i], Is.EqualTo(vec1[i]).Within(1e-12));
            }  
        }

        // UMul
        compMatrix!.UMul(compB, comp1);
        matrix.UMul(b, vec1);
        comp1.ToHost();
        {
            using var accessor = comp1.MapHost(MapFlags.Read);
            for (int i = 0; i < b.Length; i++)
            {
                Assert.That(accessor[i], Is.EqualTo(vec1[i]).Within(1e-12));
            } 
        }

        // InvUMul
        compB.CopyDeviceTo(comp1);
        b.CopyTo(vec1);

        compMatrix.InvUMul(comp1);
        matrix.InvUMul(vec1);
        comp1.ToHost();
        {
            using var accessor = comp1.MapHost(MapFlags.Read);
            for (int i = 0; i < b.Length; i++)
            {
                Assert.That(accessor[i], Is.EqualTo(vec1[i]).Within(1e-12));
            }  
        }
    }
}