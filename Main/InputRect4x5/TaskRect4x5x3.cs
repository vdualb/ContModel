#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using MathShards.Fem.Common;

public class TaskRect4x5x3: ITaskFuncs
{
    public string Description => "Прямоугольник 4на5 x^3+y^3";

    public Real Answer(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => (Real) x*x*x + y*y*y,
            _ => throw new ArgumentException("Неверный номер подобласти"),
        };
    }

    public Real Beta(int bcNum)
    {
        return bcNum switch
        {
            0 => (Real)0.5,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real F(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => (Real)(x*x*x + y*y*y - 3*(x + y)),
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Gamma(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => 1,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Lambda(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => (Real)0.5,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Theta(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => (Real)(3.0/2.0*x*x),
            1 => -(Real)(3.0/2.0*x*x),
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real uBeta(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => 3*y*y + Answer(0, x, y),
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real Ug(int bcNum, Real x, Real y)
    {
        return Answer(0, x, y);
    }
}
