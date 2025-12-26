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

using System.Text.Json;
using System.Text.Unicode;
using System.Diagnostics;

using MathShards.SlaeBuilder.Fem;
using MathShards.Matrices.Types;
using MathShards.TelmaCore;
using static MathShards.FiniteElements.Rectangle.Lagrange.BiLinear;
using MathShards.Fem.Common;
using MathShards.Mesh.RectMesh;

public class ProblemLine {
    ITaskFuncs _funcs;
    SolverParams _problemParams;
    public SolverParams ProblemParams { get => _problemParams; }
    RefineParams _refineParams;
    public FemRectMesh Mesh { get => _mesh; }
    FemRectMesh _mesh;
    MeshAxes _meshAxes;
    Subdomain[] _subDomains;
    BoundaryCondition[] _boundaryConditions;

    // TODO: подумать как сделать так, чтобы Problem
    // хранил информацию о используемом типе матрицы.
    // Сейчас настоящий тип пояляется только в момент
    // построения матрицы, в методе Build<T>().
    public IMatrix matrix;
    public Real[] b;
    public GlobalMatrixImplType buildType = GlobalMatrixImplType.Host;

    // folder - директория с условиями задачи
    public ProblemLine(
        ITaskFuncs taskFunctions,
        string taskFolder,
        GlobalMatrixImplType buildType = GlobalMatrixImplType.Host
    ) {
        _funcs = taskFunctions;

        var options = new JsonSerializerOptions()
        {
            Encoder = System.Text.Encodings.Web.JavaScriptEncoder.Create(UnicodeRanges.All),
            TypeInfoResolver = SolverParamsSourceGenerationContext.Default
        };
        
        var json = File.ReadAllText("ProblemParams.json");
        _problemParams = JsonSerializer.Deserialize<SolverParams>(json, options)!;

        options = new JsonSerializerOptions()
        {
            Encoder = System.Text.Encodings.Web.JavaScriptEncoder.Create(UnicodeRanges.All),
            TypeInfoResolver = RefineParamsSourceGenerationContext.Default
        };

        json = File.ReadAllText(Path.Combine(taskFolder, "RefineParams.json"));
        _refineParams = JsonSerializer.Deserialize<RefineParams>(json, options)!;
        // if (_refineParams.XSplitCount.Sum() <= 1)
        // {
        //     throw new Exception("Диагональный формат хранения матрицы накладывает ограничение на " +
        //         "минимальный размер сетки по оси x");
        // }

        _meshAxes = ReadMesh(taskFolder);
        _subDomains = ReadSubDomains(taskFolder);
        _boundaryConditions = ReadConditions(taskFolder);

        _mesh = new FemRectMesh(
            _meshAxes.xAxis,
            _meshAxes.yAxis,
            _subDomains,
            _boundaryConditions
        );

        _mesh.Refine(_refineParams);
        this.buildType = buildType;
    }
    
    public void Build<T>()
    where T: IFemSlaeBuilder
    {
        var sw = Stopwatch.StartNew();
        var slaeBuilder = T.Construct(_mesh, _funcs);
        slaeBuilder.GlobalMatrixImpl = buildType;
        (matrix, b) = slaeBuilder.Build();
        Trace.WriteLine($"ProblemLine.Build total: {sw.ElapsedMilliseconds}ms");
    }
    
    public void MeshRefine(RefineParams refineParams)
    {
        _mesh.Refine(refineParams);
    }

    public void MeshDouble()
    {
        _mesh.RefineDiv2();
    }

    public Real AnswerAt (Real x, Real y)
    {
        var num = _mesh.GetSubdomNumAtPoint(x, y);
        if (num.HasValue)
        {
            return _funcs.Answer(num.Value, x, y);
        }
        else
        {
            return 0;
        }
    }

    public Real ResultAt(Span<Real> q, Real x, Real y)
    {
        var X = _mesh.X;
        var Y = _mesh.Y;
        Real result = 0;

        var (xi, yi) = _mesh.GetElCoordsAtPoint(x, y);

        Real hx = X[xi + 1] - X[xi];
        Real hy = Y[yi + 1] - Y[yi];

        var subdom = _mesh.GetSubdomNumAtElCoords(xi, yi);
        if (subdom.HasValue)
        {
            Span<int> m = stackalloc int[4];
            m[0] = yi * X.Length + xi;
            m[1] = m[0] + 1;
            m[2] = (yi + 1) * X.Length + xi;
            m[3] = m[2] + 1;

            result =
                (
                      q[m[0]] * (X[xi + 1] - x) * (Y[yi + 1] - y)
                    + q[m[1]] * (x - X[xi])     * (Y[yi + 1] - y)
                    + q[m[2]] * (X[xi + 1] - x) * (y - Y[yi])
                    + q[m[3]] * (x - X[xi])     * (y - Y[yi])
                ) / hx / hy;
        }
        return result;
    }
    
    public PairReal GradAt(Span<Real> q, Real x, Real y)
    {
        var X = _mesh.X;
        var Y = _mesh.Y;
        PairReal result = new(0, 0);

        var (xi, yi) = _mesh.GetElCoordsAtPoint(x, y);

        Real hx = X[xi + 1] - X[xi];
        Real hy = Y[yi + 1] - Y[yi];

        var subdom = _mesh.GetSubdomNumAtElCoords(xi, yi);
        if (subdom.HasValue)
        {
            Span<int> m = stackalloc int[4];
            m[0] = yi * X.Length + xi;
            m[1] = m[0] + 1;
            m[2] = (yi + 1) * X.Length + xi;
            m[3] = m[2] + 1;

            var p01 = new PairReal(
                (x - X[xi]) / hx,
                (y - Y[yi]) / hy
            );

            result = new (
                q[m[0]] * BasisGrad[0, 0](p01)
                + q[m[1]] * BasisGrad[1, 0](p01)
                + q[m[2]] * BasisGrad[2, 0](p01)
                + q[m[3]] * BasisGrad[3, 0](p01)
                ,
                q[m[0]] * BasisGrad[0, 1](p01)
                + q[m[1]] * BasisGrad[1, 1](p01)
                + q[m[2]] * BasisGrad[2, 1](p01)
                + q[m[3]] * BasisGrad[3, 1](p01)
            );
            
            result = new (
                result.X / hx,
                result.Y / hy
            );
        }
        return result;
    }
    
    
    public Real ResultAtNew(Span<Real> q, Real x, Real y)
    {
        var X = _mesh.X;
        var Y = _mesh.Y;
        Real result = 0;

        var (xi, yi) = _mesh.GetElCoordsAtPoint(x, y);

        Real hx = X[xi + 1] - X[xi];
        Real hy = Y[yi + 1] - Y[yi];

        var subdom = _mesh.GetSubdomNumAtElCoords(xi, yi);
        if (subdom.HasValue)
        {
            Span<int> m = stackalloc int[4];
            m[0] = yi * X.Length + xi;
            m[1] = m[0] + 1;
            m[2] = (yi + 1) * X.Length + xi;
            m[3] = m[2] + 1;
            
            var p01 = new PairReal(
                (x - X[xi]) / hx,
                (y - Y[yi]) / hy
            );

            result =
                q[m[0]] * Basis[0](p01) 
                + q[m[1]] * Basis[1](p01)
                + q[m[2]] * Basis[2](p01)
                + q[m[3]] * Basis[3](p01);
        }
        return result;
    }

    // сохранить узлы текущего разбиения как узлы наблюдения
    // void MonitorNodesFix()
    // {
    //     // в результате выполнения этой функции у мониторных и текущих
    //     // узлов будет соотношение 1к1
    //     XMonitor = Enumerable.Range(0, slae.mesh.X.Count).ToArray();
    //     YMonitor = Enumerable.Range(0, slae.mesh.Y.Count).ToArray();
    // }

    /* Взятие нормы погрешности в пространстве Лебега 2.
        Интеграл считается методом прямоугольников  */
    public Real Lebeg2Err (Span<Real> q)
    {
        var mesh = _mesh;
        Real sum = 0;
        for (int yi = 0; yi < mesh.Y.Length - 1; yi++)
        {
            for (int xi = 0; xi < mesh.X.Length - 1; xi++)
            {
                Real x0 = mesh.X[xi];
                Real y0 = mesh.Y[yi];
                Real hx = mesh.X[xi + 1] - x0;
                Real hy = mesh.Y[yi + 1] - y0;
                var subdom = _mesh.GetSubdomNumAtElCoords(xi, yi);

                if (subdom.HasValue)
                {
                    Real u = ResultAt(q, (Real)(x0 + hx / 2d), (Real)(y0 + hy / 2d));
                    Real u_true = AnswerAt((Real)(x0 + hx / 2d), (Real)(y0 + hy / 2d));
                    Real func = u_true - u;
                    sum += hx * hy * func * func;
                }
            }
        }
        return (Real)Math.Sqrt(sum);
    }

    // public void Serialize() => matrix.Serialize();

    public (Real[] ans, int iters, Real rr) SolveOCL<T> ()
    where T: SparkAlgos.SlaeSolver.ISlaeSolver
    {
        var sw = Stopwatch.StartNew();
        var x0 = Enumerable.Repeat((Real)0, b.Length).ToArray();
        var solver = T.Construct(
            _problemParams.maxIter,
            _problemParams.eps
        );
        solver.AllocateTemps(x0.Length);
        Trace.WriteLine($"Solver prepare: {sw.ElapsedMilliseconds}ms");
        
        sw.Restart();
        
        var cm = matrix.GetComputeMatrix();
        Trace.WriteLine($"Matrix Host->Device: {sw.ElapsedMilliseconds}ms");
        
        sw.Restart();
        var (rr, iter) = solver.Solve(cm, b, x0);
        Trace.WriteLine($"Solver {typeof(T).FullName}: {sw.ElapsedMilliseconds}ms");

        return (x0, iter, rr);
    }
    
    public (Real[] ans, int iters, Real rr) SolveHost<T> ()
    where T: MathShards.SlaeSolver.ISlaeSolver
    {
        var sw = Stopwatch.StartNew();
        var x0 = Enumerable.Repeat((Real)0, b.Length).ToArray();
        var solver = T.Construct(
            _problemParams.maxIter,
            _problemParams.eps
        );
        solver.AllocateTemps(x0.Length);
        Trace.WriteLine($"Solver prepare: {sw.ElapsedMilliseconds}ms");
        
        sw.Restart();
        
        var (rr, iter) = solver.Solve(matrix, b, x0);
        Trace.WriteLine($"Solver {typeof(T).FullName}: {sw.ElapsedMilliseconds}ms");

        return (x0, iter, rr);
    }

    static MeshAxes ReadMesh(string taskFolder)
    {
        var file = new StreamReader(Path.Combine(taskFolder, "Mesh.txt"));

        MeshAxes res;
        res.xAxis = file.ReadLine()!.Split().Select(Real.Parse).ToArray();
        res.yAxis = file.ReadLine()!.Split().Select(Real.Parse).ToArray();

        return res;
    }
    
    static Subdomain[] ReadSubDomains(string taskFolder)
    {
        var file = new StreamReader(Path.Combine(taskFolder, "SubDomains.txt"));

        var domains_num = uint.Parse(file.ReadLine()!.Trim());
        var res = new Subdomain[domains_num];
        for (int i = 0; i < domains_num; i++)
        {
            var parts = file
                .ReadLine()!
                .Trim()
                .Split()
                .Select(int.Parse)
                .ToArray();
            res[i] = new Subdomain
            {
                Num = parts[0] - 1,
                X1 = parts[1] - 1,
                X2 = parts[2] - 1,
                Y1 = parts[3] - 1,
                Y2 = parts[4] - 1
            };
        }

        return res;
    }

    static BoundaryCondition[] ReadConditions(string taskFolder)
    {
        var options = new JsonSerializerOptions()
        {
            Encoder = System.Text.Encodings.Web.JavaScriptEncoder.Create(UnicodeRanges.All),
            TypeInfoResolver = BoundaryConditionsFileSourceGenerationContext.Default
        };
        
        var json = File.ReadAllText(Path.Combine(taskFolder, "BoundaryConditions.json"));
        var file = JsonSerializer.Deserialize<BoundaryConditionsFile>(json, options)!;
        
        for (int i = 0; i < file.BoundaryConditions.Length; i++)
        {
            file.BoundaryConditions[i].Num -= 1;
            file.BoundaryConditions[i].X1 -= 1;
            file.BoundaryConditions[i].X2 -= 1;
            file.BoundaryConditions[i].Y1 -= 1;
            file.BoundaryConditions[i].Y2 -= 1;
        }
        
        return file.BoundaryConditions;
    }
}
