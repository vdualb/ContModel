#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using System.Globalization;
using System.Diagnostics;

using SparkCL;
using MathShards.Matrices;
using Silk.NET.OpenCL;

namespace Main.Tests.HalvesTests;

[TestFixture]
class DiagHalvesTest
{
    const string MAIN_SRC_DIR = "../../../../Main/";
    
    [SetUp]
    public void Setup()
    {
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
        var unixMs = new DateTimeOffset(DateTime.Now).ToUnixTimeMilliseconds();
        Directory.CreateDirectory(MAIN_SRC_DIR + "measurements");
        var measures = new StreamWriter(MAIN_SRC_DIR + "measurements/" + unixMs + ".txt");
        Trace.Listeners.Add(new TextWriterTraceListener(measures));
        Trace.AutoFlush = true;
        
        Core.Init();
    }
    
    [TearDown]
    public void Teardown()
    {
        Trace.Listeners.Clear();
        
        Core.Deinit();
    }

    static (Diag9Matrix matrix, Real[] b) SomeSlae() {
        Real[] diag = [0, 0, 0, 0, 0];
        Real[] identity = [1, 1, 1, 1, 1];
        var matrix = new Diag9Matrix {
            Gap = 1,
            Di = [..identity],
            Ld3 = [..diag],
            Ld2 = [..diag],
            Ld1 = [..diag],
            Ld0 = [..diag],
            Rd3 = [..diag],
            Rd2 = [..diag],
            Rd1 = [..diag],
            Rd0 = [..diag],
        };
        matrix.Rd0[0] += 0.5f;
        matrix.Ld0[0] += 0.5f;
        Real[] b = [..identity];
        b[0] += 0.5f;        
        b[1] += 0.5f;

        return (matrix, b);
    }


    [Test]
    public static void Host()
    {
        var (matrix, b) = SomeSlae();
        Common.HalfMultiplies(matrix, b);
    }

    [Test]
    public static void OpenCL()
    {
        var (matrix, b) = SomeSlae();
        Common.HalfMultipliesOpenCL(matrix, b);
    }
}
