#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using MathShards.Fem.Common;

public class TaskRect4x5XY1: ITaskFuncs
{
    public string Description => "Прямоугольник 4на5 x+y";

    public Real Answer(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => (Real) (x + y),
            _ => throw new ArgumentException("Неверный номер подобласти"),
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

    public Real Gamma(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => 1,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real F(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => (Real)(x + y),
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Ug(int bcNum, Real x, Real y)
    {
        return Answer(0, x, y);
    }

    public Real Theta(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => 1 * Lambda(0, x, y),
            1 => -1 * Lambda(0, x, y),
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }
    
    public Real Beta(int bcNum)
    {
        return bcNum switch
        {
            0 => 0.5f,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real uBeta(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => 1 + Answer(0, x, y),
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }
}
