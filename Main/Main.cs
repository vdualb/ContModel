/*
MathShards
Copyright (C) 2025 Afonin Anton

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using System.Globalization;
using System.Diagnostics;

using SparkCL;
using MathShards.SlaeBuilder.Fem;
using MathShards.SlaeSolver;
using MathShards.CoordSystem.Dim2;
using SparkAlgos.SlaeSolver;
using MathShards.TelmaCore;
using MathShards.Mesh.RectMesh;
using MathShards.SlaeBuilder.Spline;

class Program
{
    const string SRC_DIR = "../../../";
    static void Main(string[] args)
    {
        Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
        var unixMs = new DateTimeOffset(DateTime.Now).ToUnixTimeMilliseconds();
        Directory.CreateDirectory(SRC_DIR + "measurements");
        Console.WriteLine( "Writing trace to " + Path.GetFullPath(SRC_DIR + "measurements") );
        var measures = new StreamWriter(SRC_DIR + "measurements/" + unixMs + ".txt");
        Trace.Listeners.Add(new TextWriterTraceListener(measures));
        Trace.AutoFlush = true;
        
        var sw = Stopwatch.StartNew();
        Core.Init();
        Trace.WriteLine($"SparkCL Init: {sw.ElapsedMilliseconds}ms");
        Trace.WriteLine($"Calculation type: {typeof(Real)}");

        BenchEisenstat();
        
        return;
        IndentityTest();
        HalfMultiplies();
        Spline();
        ReverseSigma();
        // ReverseI();
        ElectroMany();

        Iterate();

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
    
    static void BenchEisenstat() {
        const int REPEAT_COUNT = 5;
        const int REFINE_COUNT = 1;
        
        var task = new TaskRect4x5XY1();
        var prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");
        prob.MeshRefine(new()
        {
            XSplitCount = [1024],
            YSplitCount = [1024],
            XStretchRatio = [1],
            YStretchRatio = [1],
        });

        for (int r = 0; r < REFINE_COUNT; r++)
        {
            Console.WriteLine($"n = {prob.Mesh.nodesCount}");
            Trace.WriteLine($"n = {prob.Mesh.nodesCount}");
            prob.Build<DiagSlaeBuilder<XY>>();
            
            for (int i = 0; i < REPEAT_COUNT; i++)
            {
                SolveEisenstatHost<CgmEisenstatHost>(prob);
                // SolveEisenstatHost<CgmHost>(prob);
                // SolveOCL<CgmEisenstatOCL>(prob);
                // SolveEisenstatHost<CgmEisenstatSimpleHost>(prob);
                // SolveOCL<CgmOCL>(prob);
            }

            prob.Build<MsrSlaeBuilder>();
            
            for (int i = 0; i < REPEAT_COUNT; i++)
            {
                SolveEisenstatHost<CgmEisenstatHost>(prob);
                // SolveEisenstatHost<CgmHost>(prob);
                // SolveOCL<CgmEisenstatOCL>(prob);
                // SolveEisenstatHost<CgmEisenstatSimpleHost>(prob);
                // SolveOCL<CgmOCL>(prob);
            }

            prob.MeshDouble();
        }   
    }
    
    static void IndentityTest()
    {
        Real[] diag = [0, 0, 0, 0, 0];
        Real[] identity = [1, 1, 1, 1, 1];
        var matrix = new MathShards.Matrices.Diag9Matrix {
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
        matrix.Rd0[0] += (Real)0.5;
        matrix.Ld0[0] += (Real)0.5;
        Real[] b = [..identity];
        b[0] += (Real)0.5;        
        b[1] += (Real)0.5;        
        Real[] x = [..identity];
        x[1] += 1;
        x[2] += 1;
        x[3] += 1;
        
        var solver = CgmEisenstatHost.Construct(
            (int)1e+7,
            (Real)1e-13
        );
        solver.AllocateTemps(x.Length);
        
        var (rr, iters) = solver.Solve(matrix, b, x);
        
        Console.WriteLine("Answer: " + string.Join(", ", x));
    }
    
    static void SolveEisenstatHost<T>(ProblemLine prob)
    where T: MathShards.SlaeSolver.ISlaeSolver
    {
        Trace.WriteLine("Host solver");
        Trace.Indent();
        var sw = Stopwatch.StartNew();

        // prob.SolveHost
            var sw0 = Stopwatch.StartNew();
            var x0 = Enumerable.Repeat((Real)0, prob.b.Length).ToArray();
            var solver = T.Construct(
                prob.ProblemParams.maxIter,
                prob.ProblemParams.eps
            );
            solver.AllocateTemps(x0.Length);
            Trace.WriteLine($"Solver prepare: {sw0.ElapsedMilliseconds}ms");
            
            sw0.Restart();
            
            var (rr, iters) = solver.Solve(prob.matrix, prob.b, x0);
            var ans = x0;
            Trace.WriteLine($"Solver {typeof(T).FullName}: {sw0.ElapsedMilliseconds}ms");
        //
        
        sw.Stop();
        var err = prob.Lebeg2Err(ans);
        // Console.WriteLine($"(err {err}) (iters {iters}) (discrep: {rr})");
        Console.WriteLine($"{err} {iters} {rr}");
        Trace.WriteLine($"Solver total: {sw.ElapsedMilliseconds}мс");
        Trace.Unindent();
    }
    
    static void HalfMultiplies() {
        var task = new TaskRect4x5x3();
        var prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");

        prob.Build<DiagSlaeBuilder<XY>>();

        var (ans, iters, rr) = prob.SolveOCL<BicgStabOCL>();

        var matrix = prob.matrix as MathShards.Matrices.Diag9Matrix;

        var vec = new Real[ans.Length];
        var vec2 = new Real[ans.Length];
        ans.AsSpan().CopyTo(vec);

        matrix!.InvUMul(vec);
        matrix.UMul(vec, vec2);

        var err = vec2
            .Zip(ans)
            .Sum(a =>
            {
                return Real.Abs(a.First - a.Second);
            });
        
        err /= ans.Length;
        
        Console.WriteLine($"Error: {err}");
    }

    static (LineMesh, Real[]) CalcDiff(ProblemLine prob, Real[] q)
    {
        var basisGrad = MathShards.FiniteElements.Line.Lagrange.Linear.BasisGrad;
        
        var diffs = new Real[prob.Mesh.X.Length - 1];
        var xLine = new Real[prob.Mesh.X.Length - 1];
        for (int i = 0; i < diffs.Length; i++)
        {
            var p1 = prob.Mesh.X[i+1];
            var p0 = prob.Mesh.X[i];
            var h = p1 - p0;
            xLine[i] = (p1 + p0) / 2;
            diffs[i] = basisGrad[0](0) * q[i] + basisGrad[1](1) * q[i+1];
            diffs[i] /= h;
        }

        return (new LineMesh(xLine), diffs);
    }
    
    static Real SlopeSpline(LineMesh mesh, Real[] splineQ, Real x) {
        int c = 0;
        for (; c < mesh.X.Length; c++) {
            if (mesh.X[c+1] > x) {
                break;
            }
        }

        var x0 = mesh.X[c];
        var x1 = mesh.X[c+1];
        Real res = 0;
        var bas = MathShards.FiniteElements.Line.Hermit.Cubic.BasisConverted;
        for (int i = 0; i < 4; i++) {
            res += bas(i, x0, x1, x) * splineQ[2*c + i];
        }

        return res;
    }
    
    static Real SlopeReal(Real x) {
        return 3*x*x;
    }
    
    static Real SlopeAnsInterpol(Real[] mesh, Real[] vals, Real x) {
        int c = 0;
        for (; c < mesh.Length; c++) {
            if (mesh[c+1] > x) {
                break;
            }
        }
        
        var x0 = mesh[c];
        var x1 = mesh[c+1];
        var h = x1 - x0;
        var x01 = (x - x0) / h;
        Real res = (vals[c]*(1 - x01) + vals[c+1]*x01);

        return res;
    }
    
    static void Spline()
    {
        var task = new TaskRect4x5x3();
        var prob = new ProblemLine(task, SRC_DIR + "InputRect4x5");

        prob.Build<DiagSlaeBuilder<XY>>();

        var vals = SolveOCL<BicgStabOCL>(prob);

        // var meshSlice = new LineMesh((Real[])prob.Mesh.X.Clone());

        var (meshSlice, diffs) = CalcDiff(prob, vals);
        Console.WriteLine("Исходные значения производной");
        Console.WriteLine(string.Join(", ", meshSlice.X));
        Console.WriteLine(string.Join(", ", diffs));


        var splineProb = new ProblemSpline1D(SRC_DIR + "InputSpline", diffs, meshSlice);
        splineProb.Build<MsrSlaeBuilderHermit1D>();
        var res = SplineSolveOCL<BicgStabOCL>(splineProb);

        var splineMesh = splineProb._mesh;

        var meshOld = (Real[])meshSlice.X.Clone();
        meshSlice.RefineDiv2();
        meshSlice.RefineDiv2();
        
        Real errAns = 0;
        Real errSpline = 0;
        for (int i = 0; i < meshSlice.X.Length - 1; i++) {
            var x = meshSlice.X[i];
            var real = SlopeReal(x);
            
            var v = SlopeAnsInterpol(meshOld, diffs, x);
            errAns += Math.Abs(real - v);
            
            v = SlopeSpline(splineMesh, res, x);
            errSpline += Math.Abs(real - v);
        }

        errAns /= meshSlice.X.Length - 1;
        errSpline /= meshSlice.X.Length - 1;

        Console.WriteLine($"(errAns {errAns}) (errSpline {errSpline})");
        
        Console.WriteLine(
            string.Join(
                ",\n",
                splineProb._mesh.X.Select(
                    val => $"{val}"
                )
            )
        );

        Console.WriteLine(
            string.Join(
                ",\n",
                res.Select(
                    val => $"{val}"
                )
            )
        );
    }
    
    static Real[] SplineSolveOCL<T>(ProblemSpline prob)
    where T: SparkAlgos.SlaeSolver.ISlaeSolver
    {
        var (ans, iters, rr) = prob.SolveOCL<T>();
        
        Console.WriteLine($"(iters {iters}) (discrep: {rr})");

        return ans;
    }
    
    static Real[] SplineSolveOCL<T>(ProblemSpline1D prob)
    where T: SparkAlgos.SlaeSolver.ISlaeSolver
    {
        var (ans, iters, rr) = prob.SolveOCL<T>();
        
        Console.WriteLine($"(iters {iters}) (discrep: {rr})");

        return ans;
    }
    
    static void ReverseI()
    {
        // u -> e
        Real[] weights = [1, (Real)0.5, 1, 1];
        // измеренные напряжённости, с какой-то погрешностью
        PairReal[] e_measured = [
            new((Real)3.710883057100038e-07, -(Real)4.652739107054252e-08),
            new((Real)1.046909531845134e-07, -(Real)6.388825590065409e-09),
            new((Real)4.878075851331915e-08, -(Real)1.811845815180536e-09),
            new((Real)2.858623606945764e-08, -(Real)7.38675306794607e-10)
        ];
        // начальное значение проводимости
        Real i_guide = 2; // 1
        // истинные напряжённости
        PairReal[] e_true = [
            new((Real)3.710883057100038e-07, -(Real)4.652739107054252e-08),
            new((Real)1.046909531845134e-07, -(Real)6.388825590065409e-09),
            new((Real)4.878075851331915e-08, -(Real)1.811845815180536e-09),
            new((Real)2.858623606945764e-08, -(Real)7.38675306794607e-10)
        ];
        e_measured[1] *= (Real)1.05;
        Real alpha = 0;
        
        Real j0 = 0;
        var e0 = ElectroOnce();
        for (int i = 0; i < 4; i++)
        {
            var w = (e_measured[i] - e0[i]*i_guide) * weights[i] / e_measured[i].Norm;
            j0 += w*w;
        }
        
        Console.WriteLine($"j0: {j0}");
        
        Real beta = 1;
        Real i0 = i_guide;

        int iter = 0;
        while (true)
        {
            // дифференциал
            var diffU = new PairReal[4];
            for (int i = 0; i < 4; i++)
            {  
                diffU[i] = -e0[i];
            }
            
            // A
            Real a = alpha;
            for (int i = 0; i < 4; i++)
            {
                var w = weights[i] / e_measured[i].Norm;
                a += w*w * diffU[i]*diffU[i];
            }
            
            // f
            Real f = -alpha * (i0 - i_guide);
            for (int i = 0; i < 4; i++)
            {
                var w = weights[i] / e_measured[i].Norm;
                f -= w*w * (e_measured[i] - e0[i]*i0) * diffU[i];
            }
            
            // новое решение
            var du = f/a;
            i0 += beta * du;
            
            // новое значение функционала
            Real j1 = 0;
            for (int i = 0; i < 4; i++)
            {
                var w = (e_measured[i] - e0[i]*i0) * weights[i] / e_measured[i].Norm;
                j1 += w*w;
            }

            iter++;
            if (j0 <= j1) {
                beta /= (Real)2.0;
                if (beta < 1.0/4.0) {
                    break;
                }
            }
            
            Console.WriteLine($"I: {i0}");
            Console.WriteLine($"j: {j1}");
            
            j0 = j1;
        }
        
        Console.WriteLine($"Iterations: {iter}");
    }
    
    static void ReverseSigma()
    {
        // u -> e
        Real[] weights = [1, (Real)0.5, 1, 1];
        // измеренные напряжённости, с какой-то погрешностью
        PairReal[] e_measured = [
            new((Real)3.710883057100038e-07, -(Real)4.652739107054252e-08),
            new((Real)1.046909531845134e-07, -(Real)6.388825590065409e-09),
            new((Real)4.878075851331915e-08, -(Real)1.811845815180536e-09),
            new((Real)2.858623606945764e-08, -(Real)7.38675306794607e-10)
        ];
        // начальное значение проводимости
        Real u_guide = (Real)6.5; // 6
        // истинные напряжённости
        PairReal[] e_true = [
            new((Real)3.710883057100038e-07, -(Real)4.652739107054252e-08),
            new((Real)1.046909531845134e-07, -(Real)6.388825590065409e-09),
            new((Real)4.878075851331915e-08, -(Real)1.811845815180536e-09),
            new((Real)2.858623606945764e-08, -(Real)7.38675306794607e-10)
        ];
        e_measured[1] *= (Real)1.05;
        // e_measured = e_true;
        Real alpha = 0;
        
        Real j0 = 0;
        var e0 = ElectroOnce(u_guide);
        for (int i = 0; i < 4; i++)
        {
            var w = (e_measured[i] - e0[i]) * weights[i] / e_measured[i].Norm;
            j0 += w*w;
        }
        
        Console.WriteLine($"j0: {j0}");
        
        Real beta = 1;
        Real u0 = u_guide;

        int iter = 0;
        while (true)
        {
            // дифференциал
            var diffU = new PairReal[4];
            var e1 = ElectroOnce((Real)1.05*u0);
            for (int i = 0; i < 4; i++)
            {  
                diffU[i] = (e0[i]-e1[i])/((Real)0.05*u0);
            }
            
            // A
            Real a = alpha;
            for (int i = 0; i < 4; i++)
            {
                var w = weights[i] / e_measured[i].Norm;
                a += w*w * diffU[i]*diffU[i];
            }
            
            // f
            Real f = -alpha * (u0 - u_guide);
            for (int i = 0; i < 4; i++)
            {
                var w = weights[i] / e_measured[i].Norm;
                f -= w*w * (e_measured[i] - e0[i]) * diffU[i];
            }
            
            // новое решение
            var du = f/a;
            u0 += beta * du;
            
            // новое значение функционала
            Real j1 = 0;
            e0 = ElectroOnce(u0);
            for (int i = 0; i < 4; i++)
            {
                var w = (e_measured[i] - e0[i]) * weights[i] / e_measured[i].Norm;
                j1 += w*w;
            }

            iter++;
            if (j0 <= j1) {
                beta /= (Real)2.0;
                if (beta < 1.0/4.0) {
                    break;
                }
            }
            
            Console.WriteLine($"sigma: {u0}");
            Console.WriteLine($"j: {j1}");
            
            j0 = j1;
        }
        
        Console.WriteLine($"Iterations: {iter}");
    }
    
    static PairReal[] ElectroOnce(Real? sigma = null)
    {
        var task = new TaskElectro();
        if (sigma.HasValue) {
            task.Sigma = sigma.Value;
        }
        var prob = new ProblemLine(task, SRC_DIR + "InputElectro");
        
        // prob.MeshRefine(new()
        // {
        //     XSplitCount   = [  8,   8,   8,   8, 120],
        //     XStretchRatio = [1.0, 1.0, 1.0, 1.0, 1.0],
        //     YSplitCount   = [120,   8],
        //     YStretchRatio = [1.0, 1.0],
        // });
        
        prob.buildType = GlobalMatrixImplType.Host;
        prob.Build<DiagSlaeBuilder<RZ>>();
        // ElectroSolveHost<CgmHost>(prob);
        return ElectroSolveOCL<BicgStabOCL>(prob);
    }
    
    static void ElectroMany()
    {
        const int REFINE_COUNT = 6;
        
        var task = new TaskElectro();
        var prob = new ProblemLine(task, SRC_DIR + "InputElectro");

        for (int i = 0; i < REFINE_COUNT; i++)
        {
            prob.buildType = GlobalMatrixImplType.Host;
            prob.Build<DiagSlaeBuilder<RZ>>();
            // ElectroSolveHost<CgmHost>(prob);
            var sens = ElectroSolveOCL<BicgStabOCL>(prob);

            for (int k = 0; k < sens.Length; k++)
            {
                Console.WriteLine($"Sensor {k}: {sens[k]}.");
            }

            prob.MeshDouble();
        }
    }
    
    static PairReal[] ElectroSolveOCL<T>(ProblemLine prob)
    where T: SparkAlgos.SlaeSolver.ISlaeSolver
    {
        var (ans, iters, rr) = prob.SolveOCL<T>();
        
        Console.WriteLine($"(iters {iters}) (discrep: {rr})");

        // sensors
        Real[] sensR = [250, 500, 750, 1000];
        
        var sensVals = new PairReal[4];
        var res = new Real[4];
        var resNew = new Real[4];
        for (int i = 0; i < sensR.Length; i++)
        {
            sensVals[i] = -prob.GradAt(ans, sensR[i], 0);
            res[i] = prob.ResultAt(ans, sensR[i], 0);
            resNew[i] = prob.ResultAtNew(ans, sensR[i], 0);
        }
        return sensVals;
    }
    
    static void ElectroSolveHost<T>(ProblemLine prob)
    where T: MathShards.SlaeSolver.ISlaeSolver
    {
        var (ans, iters, rr) = prob.SolveHost<T>();
        
        Console.WriteLine($"(iters {iters}) (discrep: {rr})");

        // sensors
        Real[] sensR = [250, 500, 750, 1000];

        for (int i = 0; i < sensR.Length; i++)
        {
            // var res = ans[sensIdx[i]];
            var res = prob.ResultAt(ans, sensR[i], 0);
            Console.WriteLine($"Sensor {i+1}: {res}");
        }

        Console.WriteLine();
    }
    
    static void Iterate()
    {
        const int REFINE_COUNT = 3;

        var task = new TaskRect4x5RZ2();
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
            prob.Build<DiagSlaeBuilder<RZ>>();
            SolveOCL<BicgStabOCL>(prob);

            prob.MeshDouble();
        }
    }
    
    static Real[] SolveOCL<T>(ProblemLine prob)
    where T: SparkAlgos.SlaeSolver.ISlaeSolver
    {
        Trace.WriteLine("OpenCL solver");
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
        // Console.WriteLine($"(err {err}) (iters {iters}) (discrep: {rr})");
        Console.WriteLine($"{err} {iters} {rr}");
        Trace.Unindent();
        Trace.Write($"Solver total: {sw.ElapsedMilliseconds}мс");
#if SPARKCL_COLLECT_TIME
        Trace.Write($": {kernTime}мс + {ioTime}мс");
#endif
        Trace.WriteLine("");

        return ans;
    }
}
