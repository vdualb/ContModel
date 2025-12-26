#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using MathShards.Fem.Common;

public class TaskRect4x5Exp : ITaskFuncs
{
    public string Description => "Прямоугольник 4на5 с экспонентой";

    public Real Answer(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => Real.Exp(x+y) + x,
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
            0 => (y*y-1)*Real.Exp(x+y) + x*y*y,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Gamma(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => y*y,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Lambda(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => 0.5f,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Theta(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => (Real.Exp(4 + y) + 1)/2,
            1 => -(Real.Exp(1+y) + 1)/2,
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real uBeta(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => 2*Real.Exp(x+5) + x,
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real Ug(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => Real.Exp(x+y) + x,
            1 => Real.Exp(x+y) + x,
            2 => Real.Exp(x+y) + x,
            3 => Real.Exp(x+y) + x,
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

}
