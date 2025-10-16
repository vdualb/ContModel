using Real = double;

using TelmaCore;
using static Quadrature.Gauss;

namespace FiniteElements.Line.Hermit;

public static class Cubic
{
    public static readonly Func<Real, Real>[] Basis =
    {
        a => 1 - 3*a*a + 2*a*a*a,
        a => a - 2*a*a + a*a*a,
        a => 3*a*a - 2*a*a*a,
        a => -a*a + a*a*a
    };
}
