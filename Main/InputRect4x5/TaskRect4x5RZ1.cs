#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using MathShards.Fem.Common;

public class TaskRect4x5RZ1: ITaskFuncs
{
    public string Description => "Прямоугольник 4на5 в цилиндрических координатах, полинов второй степени";

    public Real Answer(int subdom, Real r, Real z)
    {
        return subdom switch
        {
            0 => r + z,
            _ => throw new ArgumentException("Неверный номер подобласти"),
        };
    }

    public Real Beta(int bcNum)
    {
        return bcNum switch
        {
            0 => (Real)0.5,
            1 => (Real)0.5,
            2 => (Real)0.5,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real F(int subdom, Real r, Real z)
    {
        return subdom switch
        {
            0 => r + z - (Real)(1.0/2.0/r),
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Gamma(int subdom, Real r, Real z)
    {
        return subdom switch
        {
            0 => 1,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Lambda(int subdom, Real r, Real z)
    {
        return subdom switch
        {
            0 => (Real)0.5,
            _ => throw new ArgumentException("Неверный номер граничного условия"),
        };
    }

    public Real Theta(int bcNum, Real r, Real z)
    {
        return bcNum switch
        {
            0 => 1 * Lambda(0, r, z),
            1 => -1 * Lambda(0, r, z),
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real uBeta(int bcNum, Real r, Real z)
    {
        return bcNum switch
        {
            0 => 1 + Answer(0, r, z),
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real Ug(int bcNum, Real r, Real z)
    {
        return Answer(0, r, z);
    }
}
