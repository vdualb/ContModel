#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using System.Globalization;
using System.Diagnostics;

using SparkCL;
using SlaeBuilder;
using SparkAlgos.SlaeSolver;
using SlaeSolver;

class Program
{
    const string SRC_DIR = "../../../";
    static void Main(string[] args)
    {
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
        var unixMs = new DateTimeOffset(DateTime.Now).ToUnixTimeMilliseconds();
        Directory.CreateDirectory(SRC_DIR + "measurements");
        var measures = new StreamWriter(SRC_DIR + "measurements/" + unixMs + ".txt");
        Trace.Listeners.Add(new TextWriterTraceListener(measures));
        Trace.AutoFlush = true;
        
        var sw = Stopwatch.StartNew();
        Core.Init();
        Trace.WriteLine($"SparkCL Init: {sw.ElapsedMilliseconds}ms");

        ElectroMany();

        return;
        Reverse();


        var e = ElectroOnce(7);
        Console.WriteLine(
            string.Join(
                "\n",
                e.Select(
                    (val, idx) => $"{idx}: {val}."
                )
            )
        );
    }
    
    static void Reverse()
    {
        // u -> e
        Real[] w = [1, 1, 1, 1];
        Real[] e_measured = [37e-3, 18e-3, 12e-3, 8e-3];
        Real u_guide = 6; // 7
        Real[] e_true = [
            0.03757434500495688,
            0.01803733811906942,
            0.0117831743630933,
            0.008658142365752039
        ];
        e_measured = e_true;
        Real alpha = 0.01;
        
        Real j0 = 0;
        var e0 = ElectroOnce(u_guide);
        for (int i = 0; i < 4; i++)
        {
            j0 += w[i]*w[i] * (e_measured[i] - e0[i]) * (e_measured[i] - e0[i]);
        }
        
        Console.WriteLine($"j0: {j0}");
        
        Real beta = 1;
        Real u0 = u_guide;
        
        while (true)
        {
            // дифференциал
            var diffU = new Real[4];
            var e1 = ElectroOnce(1.05*u0);
            for (int i = 0; i < 4; i++)
            {  
                diffU[i] = (e1[i]-e0[i])/(0.05*u0);
            }
            
            // A
            Real a = alpha;
            for (int i = 0; i < 4; i++)
            {
                a += w[i]*w[i] * diffU[i]*diffU[i];
            }
            
            // f
            Real f = -alpha * (u0 - u_guide);
            for (int i = 0; i < 4; i++)
            {
                f -= w[i]*w[i] * (e_measured[i] - e0[i])  * diffU[i];
            }
            
            // новое решение
            var du = f/a;
            u0 += beta * du;
            
            // новое значение функционала
            Real j1 = 0;
            e0 = ElectroOnce(u0);
            for (int i = 0; i < 4; i++)
            {
                j1 += w[i]*w[i] * (e_measured[i] - e0[i]) * (e_measured[i] - e0[i]);
            }
            
            if (j0 < j1) {
                beta /= 2;
                if (beta < 1e-7) {
                    break;
                }
            }
            
            Console.WriteLine($"sigma: {u0}");
            Console.WriteLine($"j: {j1}");
            
            j0 = j1;
        }
    }
    
    static Real[] ElectroOnce(Real sigma)
    {
        var task = new TaskElectro();
        task.Sigma = sigma;
        var prob = new ProblemLine(task, SRC_DIR + "InputElectro");
        
        prob.MeshRefine(new()
        {
            XSplitCount   = [  8,   8,   8,   8, 120],
            XStretchRatio = [1.0, 1.0, 1.0, 1.0, 1.0],
            YSplitCount   = [120,   8],
            YStretchRatio = [1.0, 1.0],
        });
        
        prob.buildType = GlobalMatrixImplType.Host;
        prob.Build<DiagSlaeBuilder>();
        // ElectroSolveHost<CgmHost>(prob);
        return ElectroSolveOCL<CgmOCL>(prob);
    }
    
    static void ElectroMany()
    {
        const int REFINE_COUNT = 6;
        
        var task = new TaskElectro();
        var prob = new ProblemLine(task, SRC_DIR + "InputElectro");

        for (int i = 0; i < REFINE_COUNT; i++)
        {
            prob.buildType = GlobalMatrixImplType.Host;
            prob.Build<DiagSlaeBuilder>();
            // ElectroSolveHost<CgmHost>(prob);
            var sens = ElectroSolveOCL<CgmOCL>(prob);

            for (int k = 0; k < sens.Length; k++)
            {
                Console.WriteLine($"Sensor {k}: {sens[k]}.");
            }

            prob.MeshDouble();
        }
    }
    
    static Real[] ElectroSolveOCL<T>(ProblemLine prob)
    where T: SparkAlgos.SlaeSolver.ISlaeSolver
    {
        var (ans, iters, rr) = prob.SolveOCL<T>();
        
        Console.WriteLine($"(iters {iters}) (discrep: {rr})");

        // sensors
        int[] sensIdx = [
            prob.Mesh.GetDofAtInitNode(1, 2),
            prob.Mesh.GetDofAtInitNode(2, 2),
            prob.Mesh.GetDofAtInitNode(3, 2),
            prob.Mesh.GetDofAtInitNode(4, 2)
        ];
        
        var sensVals = new Real[4];

        for (int i = 0; i < sensIdx.Length; i++)
        {
            sensVals[i] = ans[sensIdx[i]];
        }
        return sensVals;
    }
    
    static void ElectroSolveHost<T>(ProblemLine prob)
    where T: SlaeSolver.ISlaeSolver
    {
        var (ans, iters, rr) = prob.SolveHost<T>();
        
        Console.WriteLine($"(iters {iters}) (discrep: {rr})");

        // sensors
        int[] sensIdx = [
            prob.Mesh.GetDofAtInitNode(1, 2),
            prob.Mesh.GetDofAtInitNode(2, 2),
            prob.Mesh.GetDofAtInitNode(3, 2),
            prob.Mesh.GetDofAtInitNode(4, 2)
        ];

        for (int i = 0; i < sensIdx.Length; i++)
        {
            var res = ans[sensIdx[i]];
            Console.WriteLine($"Sensor {i+1}: {res}");
        }

        Console.WriteLine();
    }
    
    static void Iterate()
    {
        const int REFINE_COUNT = 3;

        var task = new TaskRect4x5RZ();
        var prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");

        // prob.MeshRefine(new()
        // {
        //     XSplitCount = [32],
        //     YSplitCount = [32],
        //     XStretchRatio = [1],
        //     YStretchRatio = [1],
        // });

        for (int i = 0; i < REFINE_COUNT; i++)
        {
            prob.buildType = GlobalMatrixImplType.Host;
            prob.Build<DiagSlaeBuilder>();
            SolveOCL<BicgStabOCL>(prob);

            prob.MeshDouble();
        }
    }
    
    static void SolveOCL<T>(ProblemLine prob)
    where T: SparkAlgos.SlaeSolver.ISlaeSolver
    {
        Trace.WriteLine("SolveOCL");
        Trace.Indent();
#if SPARKCL_COLLECT_TIME
        Core.ResetTime();
#endif
        var sw = Stopwatch.StartNew();
        var (ans, iters, rr) = prob.SolveOCL<T>();
        sw.Stop();
        var err = prob.Lebeg2Err(ans);
#if SPARKCL_COLLECT_TIME
        var (ioTime, kernTime) = Core.MeasureTime();
        ioTime /= (ulong)1e+6;
        kernTime /= (ulong)1e+6;
#endif
        Console.WriteLine($"(err {err}) (iters {iters}) (discrep: {rr})");
        Trace.Unindent();
        Trace.Write($"Solver total: {sw.ElapsedMilliseconds}мс");
#if SPARKCL_COLLECT_TIME
        Trace.Write($": {kernTime}мс + {ioTime}мс");
#endif
        Trace.WriteLine("");
    }
}
