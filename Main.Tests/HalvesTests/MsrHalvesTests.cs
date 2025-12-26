#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using System.Globalization;
using System.Diagnostics;

using SparkCL;
using MathShards.CoordSystem.Dim2;
using MathShards.SlaeBuilder.Fem;
using SparkAlgos.SlaeSolver;
using MathShards.Matrices;
using Silk.NET.OpenCL;

namespace Main.Tests.HalvesTests;

[TestFixture]
class MsrHalvesTest
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

    static (MsrMatrix matrix, Real[] b) SomeSlae() {
        Real[] identity = [1, 1, 1, 1, 1];
        Real[] elems = [0.5f, 0.5f, 0.5f, 0.5f];

        var matrix = new MsrMatrix {
            Elems = [..elems],
            Ia = [0, 1, 2, 3, 3, 4],
            Ja = [1, 0, 3, 3],
            Di = [..identity],
        };
        Real[] b = [1.5f, 1.5f, 1, 1, 1];

        return (matrix, b);
    } 

    static (MsrMatrix matrix, Real[] b) ComplexSlae() {
        Real[] elems = [0.5f, 0.5f, 0.2f, 0.3f, 0.2f, 0.5f, 0.3f, 0.5f];
        var matrix = new MsrMatrix {
            Elems = [..elems],
            Ia = [0, 1, 4, 6, 7, 8],
            Ja = [1, 0, 2, 3, 1, 4, 1, 2],
            Di = [1.5f, 2, 2.5f, 3, 3.5f],
        };
        Real[] b = [0.5f, 1, 1.5f, 2, 3.5f];

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
