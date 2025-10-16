using Real = double;

using System.Diagnostics;

using Matrices;
using Types;
using Dim2 = FiniteElements.Rectangle.Lagrange.BiLinear;
using static Quadrature.Gauss;
using TelmaCore;

namespace SlaeBuilder;

class DiagSlaeBuilderHermit : ISlaeBuilder
{
    Diag9Matrix _matrix;
    Real[] _b = [];

    readonly RectMesh _mesh;
    public RectMesh Mesh { get => _mesh; }
    public GlobalMatrixImplType GlobalMatrixImpl { get; set; } = GlobalMatrixImplType.Host;

    readonly TaskFuncs _funcs;

    public DiagSlaeBuilderHermit(RectMesh mesh, TaskFuncs funcs)
    {
        throw new NotImplementedException();
        _mesh = mesh;
        _matrix = new Diag9Matrix();
        _funcs = funcs;
    }

    public static ISlaeBuilder Construct(RectMesh mesh, TaskFuncs funcs)
        => new DiagSlaeBuilderHermit(mesh, funcs);

    public (IMatrix, Real[]) Build()
    {
        Trace.WriteLine($"Diag Builder: {GlobalMatrixImpl}");

        Trace.Indent();
        var sw = Stopwatch.StartNew();
        GlobalMatrixInit();
        Trace.WriteLine($"Init: {sw.ElapsedMilliseconds}ms");

        sw.Restart();
        GlobalMatrixBuild();
        Trace.WriteLine($"Build: {sw.ElapsedMilliseconds}ms");

        return (_matrix, _b);
    }

    void GlobalMatrixInit()
    {
        GlobalMatrixPortraitCompose();

        _matrix.Di = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();

        _matrix.Ld0 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.Ld1 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.Ld2 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.Ld3 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();

        _matrix.Rd0 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.Rd1 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.Rd2 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
        _matrix.Rd3 = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();

        _b = Enumerable.Repeat((Real)0, _mesh.nodesCount).ToArray();
    }

    void GlobalMatrixPortraitCompose()
    {
        // в случае диагональной матрицы это просто Gap
        _matrix.Gap = _mesh.X.Length - 2;
    }

    void GlobalMatrixBuild()
    {
        switch (GlobalMatrixImpl)
        {
            case GlobalMatrixImplType.Host:
                GlobalMatrixBuildImplHost();
                break;
            default:
                throw new InvalidOperationException();
        }
    }

    void GlobalMatrixBuildImplHost()
    {
        for (int yi = 0; yi < _mesh.Y.Length - 1; yi++)
        {
            for (int xi = 0; xi < _mesh.X.Length - 1; xi++)
            {
                var subDom = _mesh.GetSubdomNumAtElCoords(xi, yi);
                if (subDom == null) continue;
                
                PairF64 p0 = new(_mesh.X[xi], _mesh.Y[yi]);
                PairF64 p1 = new(_mesh.X[xi + 1], _mesh.Y[yi + 1]);
                
                var local = Dim2.ComputeLocal(_funcs, p0, p1, subDom.Value);

                int a = yi * _mesh.X.Length + xi;

                _matrix.Di[a] += local[0, 0];
                _matrix.Ld0[a] += local[0, 1];
                _matrix.Rd0[a] += local[0, 1];
                _matrix.Ld2[a] += local[0, 2];
                _matrix.Rd2[a] += local[0, 2];
                _matrix.Ld3[a] += local[0, 3];
                _matrix.Rd3[a] += local[0, 3];

                _matrix.Di[a + 1] += local[1, 1];
                _matrix.Ld1[a + 1] += local[1, 2];
                _matrix.Rd1[a + 1] += local[1, 2];
                _matrix.Ld2[a + 1] += local[1, 3];
                _matrix.Rd2[a + 1] += local[1, 3];

                var a2 = a + _mesh.X.Length;
                _matrix.Di[a2] += local[2, 2];
                _matrix.Ld0[a2] += local[2, 3];
                _matrix.Rd0[a2] += local[2, 3];

                _matrix.Di[a2 + 1] += local[3, 3];

                var localB = Dim2.ComputeLocalB(_funcs, p0, p1, subDom.Value);

                _b[a] += localB[0];
                _b[a + 1] += localB[1];
                _b[a2] += localB[2];
                _b[a2 + 1] += localB[3];
            }
        }

        /* После сборки матрицы надо нулевые диагональные элементы заменить
            на 1 */
        for (int i = 0; i < _matrix.Di.Length; i++)
        {
            if (_matrix.Di[i] == 0)
            {
                _matrix.Di[i] = 1;
            }
        }
    }
}
