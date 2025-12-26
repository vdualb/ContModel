#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using static MathShards.Quadrature.Gauss;
using MathShards.TelmaCore;

namespace Main.Tests;

[TestFixture]
class QuadratureTests
{
    [SetUp]
    public void Setup()
    {
    }

    [Test]
    public void Gauss31DTemplateElement()
    {
        Real p0 = -1;
        Real p1 = 1;

        var func = (Real x) => x * x * x + 2 * x + 1;
        var res = Integrate1DOrder5(p0, p1, func);
        
        Assert.That(res, Is.EqualTo(2));        
    }
    
    [Test]
    public void Gauss31DNonTemplate()
    {
        Real p0 = -10;
        Real p1 = 20;

        static Real func(Real x) => (Real)(0.1 * x * x - 2 * x - 10);
        var res = Integrate1DOrder5(p0, p1, func);
        
        Assert.That(res, Is.EqualTo(-300));
    }
    
    
    [Test]
    public void Gauss32DTemplate() {
        PairReal p0 = new(-1, -1);
        PairReal p1 = new(1, 1);

        var func = (PairReal p) => p.X * p.X * p.Y * p.Y + p.X + p.Y;
        var res = Integrate2DOrder5(p0, p1, func);
        
        #if USE_DOUBLE
        Assert.That(res, Is.EqualTo(4.0 / 9.0).Within(1e-13));
        #else
        Assert.That(res, Is.EqualTo(4.0 / 9.0).Within(1e-6));
        #endif
    }
    
    [Test]
    public void Gauss32DNonTemplate() {
        PairReal p0 = new(-10, -10);
        PairReal p1 = new(10, 10);

        var func = (PairReal p) => p.X * p.X * p.Y * p.Y + p.X + p.Y;
        var res = Integrate2DOrder5(p0, p1, func);
        
        #if USE_DOUBLE
        Assert.That(res, Is.EqualTo(444444 + 4.0 / 9.0).Within(1e-7));
        #else
        Assert.That(res, Is.EqualTo(444444 + 4.0 / 9.0).Within(1));
        #endif
    }
}
