using Real = double;

using TelmaCore;
using static Quadrature.Gauss;

namespace FiniteElements.Line.Lagrange;

public static class Linear
{
    public static readonly Func<Real, Real>[] LagrangeBasis =
    {
        coord => 1 - coord,
        coord => coord
    };
}
