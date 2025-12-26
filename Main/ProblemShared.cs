/*
MathShards
Copyright (C) 2025 Afonin Anton

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

using System.Text.Json.Serialization;
using MathShards.Fem.Common;
using MathShards.Mesh.RectMesh;

[JsonSerializable(typeof(SolverParams))]
[JsonSerializable(typeof(Real))]
[JsonSerializable(typeof(int))]        
internal partial class SolverParamsSourceGenerationContext : JsonSerializerContext
{
}

public struct SolverParams
{
    public Real eps { get; set; }
    public int maxIter { get; set; }
}

[JsonSerializable(typeof(RefineParams))]
[JsonSerializable(typeof(int[]))]
[JsonSerializable(typeof(Real[]))]        
internal partial class RefineParamsSourceGenerationContext : JsonSerializerContext
{
}


[JsonSerializable(typeof(int))]
[JsonSerializable(typeof(BoundaryCondition[]))]        
[JsonSerializable(typeof(BoundaryConditionsFile))]        
internal partial class BoundaryConditionsFileSourceGenerationContext
    : JsonSerializerContext
{}

public struct BoundaryConditionsFile
{
    public BoundaryCondition[] BoundaryConditions { get; set; }
}

public struct MeshAxes
{
    public Real[] xAxis;
    public Real[] yAxis;
}

interface IElement
{
    int NumberOfDofs { get; }
    Real[,] LocalGMatrix (Real hx, Real hy, Real gamma);
    Real[,] LocalMMatrix (Real hx, Real hy, Real gamma);
}
