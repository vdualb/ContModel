using Real = double;

using TelmaCore;
using static Quadrature.Gauss;

namespace FiniteElements.Rectangle.Hermit;
using Dim1 = FiniteElements.Line.Hermit.Cubic;

public static class Cubic
{
    public static Func<PairF64, Real> Basis(int i)
    {
        int mu = 2*((i/4)%2) +i%2;
        int nu = 2*(i/8) + i/2%2;

        return pair
            => Dim1.Basis[mu](pair.X) * Dim1.Basis[nu](pair.Y);
    }
}
