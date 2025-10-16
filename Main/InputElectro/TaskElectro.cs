using Real = double;

class TaskElectro : TaskFuncs
{
    public string Description => "Двухслойная среда";

    public Real Answer(int subdom, Real x, Real y)
    {
        throw new NotSupportedException();
    }

    public Real Sigma { get; set; } = 7;
    public Real Lambda(int subdom, Real x, Real y)
    {
        return subdom switch
        {
            0 => 6,
            1 => Sigma,
            _ => throw new ArgumentException("Неверный номер области"),
        };
    }
    
    public Real Gamma(int subdom, Real x, Real y)
    {
        return 0;
    }

    public Real F(int subdom, Real x, Real y)
    {
        return 0;
    }

    public Real Ug(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => 0,
            1 => 0,
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }

    public Real Theta(int bcNum, Real x, Real y)
    {
        return bcNum switch
        {
            0 => 0,
            1 => 0,
            _ => throw new ArgumentException("Некорректный номер условия"),
        };
    }
    
    public Real Beta(int bcNum)
    {
        throw new NotSupportedException();
    }

    public Real uBeta(int bcNum, Real x, Real y)
    {
        throw new NotSupportedException();
    }
}
