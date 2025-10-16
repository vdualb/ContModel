using static Quadrature.Gauss;
using TelmaCore;

namespace Main.Tests;

public class Tests
{
    [SetUp]
    public void Setup()
    {
    }

    [Test]
    public void Gauss31DTemplateElement()
    {
        double p0 = -1;
        double p1 = 1;

        var func = (double x) => x * x * x + 2 * x + 1;
        var res = Integrate1DOrder5(p0, p1, func);
        
        Assert.That(res, Is.EqualTo(2));        
    }
    
    [Test]
    public void Gauss31DNonTemplate()
    {
        double p0 = -10;
        double p1 = 20;

        var func = (double x) => 0.1 * x * x - 2 * x - 10;
        var res = Integrate1DOrder5(p0, p1, func);
        
        Assert.That(res, Is.EqualTo(-300));
    }
    
    
    [Test]
    public void Gauss32DTemplate() {
        PairF64 p0 = new(-1, -1);
        PairF64 p1 = new(1, 1);

        var func = (PairF64 p) => p.X * p.X * p.Y * p.Y + p.X + p.Y;
        var res = Integrate2DOrder5(p0, p1, func);
        
        Console.WriteLine(res);

        Assert.That(res, Is.EqualTo(4.0 / 9.0).Within(1e-13));
    }
    
    [Test]
    public void Gauss32DNonTemplate() {
        PairF64 p0 = new(-10, -10);
        PairF64 p1 = new(10, 10);

        var func = (PairF64 p) => p.X * p.X * p.Y * p.Y + p.X + p.Y;
        var res = Integrate2DOrder5(p0, p1, func);
        
        Console.WriteLine(res);

        Assert.That(res, Is.EqualTo(444444 + 4.0 / 9.0).Within(1e-7));
    }
}
