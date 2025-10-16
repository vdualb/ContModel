using Real = double;

using TelmaCore;
using static Quadrature.Gauss;

namespace FiniteElements.Rectangle.Lagrange;

public static class BiLinear
{
    public static readonly Func<PairF64, Real>[] Basis =
    {
        vert => (1 - vert.X) * (1 - vert.Y),
        vert => vert.X       * (1 - vert.Y),
        vert => (1 - vert.X) * vert.Y,
        vert => vert.X       * vert.Y
    };

    public static readonly Func<PairF64, Real>[,] BasisGrad =
    {
        {
            vert => -(1 - vert.Y),
            vert => -(1 - vert.X),
        },
        {
            vert => (1 - vert.Y),
            vert => -vert.X,
        },
        {
            vert => -vert.Y,
            vert => (1 - vert.X),
        },
        {
            vert => vert.Y,
            vert => vert.X
        }
    };
    
    public static Real[,] ComputeLocal(TaskFuncs funcs, PairF64 p0, PairF64 p1, int subDom)
    {
        var values = new Real[4, 4];

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                var ph = p1 - p0;
                
                var funcMass = (PairF64 point) => {
                    // в координатах шаблонного базиса, [0;1]
                    var p01 = new PairF64 (
                        (point.X - p0.X) / ph.X,
                        (point.Y - p0.Y) / ph.Y
                    );
                    return funcs.Gamma(subDom, point.X, point.Y)
                        * Basis[i](p01)
                        * Basis[j](p01)
                        * point.X;
                };
                
                values[i, j] = Integrate2DOrder5(p0, p1, funcMass);
                
                var funcStiffness = (PairF64 point) => {
                    // в координатах шаблонного базиса - [0;1]
                    var p01 = new PairF64 (
                        (point.X - p0.X) / ph.X,
                        (point.Y - p0.Y) / ph.Y
                    );
                    return  funcs.Lambda(subDom, point.X, point.Y)
                        *
                        (
                            BasisGrad[i, 0](p01)
                            * BasisGrad[j, 0](p01) / ph.X / ph.X
                        +
                            BasisGrad[i, 1](p01)
                            * BasisGrad[j, 1](p01) / ph.Y / ph.Y
                        )
                        * point.X;
                };
                
                values[i, j] += Integrate2DOrder5(p0, p1, funcStiffness);
            }
        }

        return values;
    }
    
    public static Real[] ComputeLocalB(TaskFuncs funcs, PairF64 p0, PairF64 p1, int subDom)
    {
        var ph = p1 - p0;
        var res = new Real[4];
        
        for (int i = 0; i < 4; i++)
        {
            var func = (PairF64 point) =>
            {
                // в координатах шаблонного базиса - [0;1]
                var p01 = new PairF64 (
                    (point.X - p0.X) / ph.X,
                    (point.Y - p0.Y) / ph.Y
                );
                return funcs.F(subDom, point.X, point.Y)
                    * Basis[i](p01)
                    * point.X;
            };
            res[i] = Integrate2DOrder5(p0, p1, func);
        }
        
        return res;
    }
    
    // для декартовой системы координат
    public static readonly Real[,] LocalG1 = {
        { 2, -2,  1, -1},
        {-2,  2, -1,  1},
        { 1, -1,  2, -2},
        {-1,  1, -2,  2},
    };
    public static readonly Real[,] LocalG2 = {
        { 2,  1, -2, -1},
        { 1,  2, -1, -2},
        {-2, -1,  2,  1},
        {-1, -2,  1,  2},
    };
    public static readonly Real[,] LocalM = {
        {4, 2, 2, 1},
        {2, 4, 1, 2},
        {2, 1, 4, 2},
        {1, 2, 2, 4},
    };
    
    public static Real[] ComputeLocalBTempl(TaskFuncs funcs, PairF64 p0, PairF64 p1, int subDom)
    {
        var ph = p1 - p0;
        var res = new Real[4];
        
        /* правая часть */
        Real f1 = funcs.F(subDom, p0.X, p0.Y);
        Real f2 = funcs.F(subDom, p1.X, p0.Y);
        Real f3 = funcs.F(subDom, p0.X, p1.Y);
        Real f4 = funcs.F(subDom, p1.X, p1.Y);

        Real hx = p1.X - p0.X;
        Real hy = p1.Y - p0.Y;

        res[0] = hx * hy / 36 * (4 * f1 + 2 * f2 + 2 * f3 + f4);
        res[1] = hx * hy / 36 * (2 * f1 + 4 * f2 + f3 + 2 * f4);
        res[2] = hx * hy / 36 * (2 * f1 + f2 + 4 * f3 + 2 * f4);
        res[3] = hx * hy / 36 * (f1 + 2 * f2 + 2 * f3 + 4 * f4);
    
        return res;
    }

    public static Real [,] ComputeLocalTempl(TaskFuncs funcs, PairF64 p0, PairF64 p1, int subDom)
    {
        Real GetGammaAverage()
        {
            Real res = funcs.Gamma(subDom, p0.X, p0.Y)
                     + funcs.Gamma(subDom, p1.X, p0.Y)
                     + funcs.Gamma(subDom, p0.X, p1.Y)
                     + funcs.Gamma(subDom, p1.X, p1.Y);

            return res / 4;
        }

        Real GetLamdaAverage()
        {
            Real res = funcs.Lambda(subDom, p0.X, p0.Y)
                     + funcs.Lambda(subDom, p1.X, p0.Y)
                     + funcs.Lambda(subDom, p0.X, p1.Y)
                     + funcs.Lambda(subDom, p1.X, p1.Y);

            return res / 4;
        }
        
        Real hy = p1.Y - p0.Y;
        Real hx = p1.X - p0.X;
        Real l_avg = GetLamdaAverage();
        Real g_avg = GetGammaAverage();
        
        var values = new Real[4, 4];
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                values[i, j] = l_avg / 6 * (hy / hx * LocalG1[i, j] + hx / hy * LocalG2[i, j])
                    + g_avg / 36 * hx * hy * LocalM[i, j];
            }
        }

        return values;
    }
}
