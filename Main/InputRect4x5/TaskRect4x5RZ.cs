using Real = double;

class TaskRect4x5RZ : TaskFuncs
{
    public string Description => "Прямоугольник 4на5 в цилиндрических координатах";

    public Real Answer(int subdom, Real r, Real z)
    {
        return subdom switch
        {
            0 => r*r + z*z,
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
            0 => r*r + z*z - 3.0,
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
            0 => 4,
            1 => -1,
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real uBeta(int bcNum, Real r, Real z)
    {
        return bcNum switch
        {
            0 => 2*z + Answer(0, r, z),
            1 => -2*r + Answer(0, r, z),
            2 => 2*r + Answer(0, r, z),
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real Ug(int bcNum, Real r, Real z)
    {
        return Answer(0, r, z);
    }
}
