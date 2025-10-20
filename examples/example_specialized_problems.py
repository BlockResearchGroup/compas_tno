"""
Example: Using Specialized Problem Classes
===========================================

This example demonstrates the new specialized Problem classes in compas_tno.
"""

from compas_tna.diagrams import FormDiagram
from compas_tno.problems import FixedProblem
from compas_tno.problems import FixedSymmetricProblem
from compas_tno.problems import Problem
from compas_tno.problems import SymmetricProblem

# =============================================================================
# Example 1: Base Problem (No constraints)
# =============================================================================

print("=" * 70)
print("Example 1: Base Problem - All edges are independent")
print("=" * 70)

x_span = (0, 10.0)
y_span = (0, 10.0)
n = 5

form = FormDiagram.create_cross(x_span=x_span, y_span=y_span, n=n)
problem = Problem.from_formdiagram(form)

print(f"Vertices: {problem.n}")
print(f"Edges: {problem.m}")
print(f"Independent variables: {problem.k}")
print(f"All edges are independent: {problem.k == problem.m}")
print()

# =============================================================================
# Example 2: Fixed Problem (Diagram fixed in plan)
# =============================================================================

print("=" * 70)
print("Example 2: FixedProblem - Independent edges computed")
print("=" * 70)

form = FormDiagram.create_cross(x_span=x_span, y_span=y_span, n=n)
problem = FixedProblem.from_formdiagram(form, method="QR", printout=False)

print(f"Vertices: {problem.n}")
print(f"Edges: {problem.m}")
print(f"Independent variables: {problem.k}")
print(f"Dependent variables: {len(problem.dep)}")
print(f"Reduction: {100 * (1 - problem.k / problem.m):.1f}%")
print()

# =============================================================================
# Example 3: Symmetric Problem (Exploits symmetry)
# =============================================================================

print("=" * 70)
print("Example 3: SymmetricProblem - Symmetry reduces variables")
print("=" * 70)

center = (5.0, 5.0)
radius = 5.0
n_hoops = 6
n_parallels = 12
r_oculus = 0.5

form = FormDiagram.create_circular_radial_spaced(
    center=center,
    radius=radius,
    n_hoops=n_hoops,
    n_parallels=n_parallels,
    r_oculus=r_oculus,
)

# Remove boundary edges
for edge in form.edges_on_boundary():
    form.edge_attribute(edge, "_is_edge", False)

problem = SymmetricProblem.from_formdiagram(form, center=center, printout=False)

print(f"Vertices: {problem.n}")
print(f"Edges: {problem.m}")
print(f"Independent variables: {problem.k}")
print(f"Reduction through symmetry: {100 * (1 - problem.k / problem.m):.1f}%")
print()

# =============================================================================
# Example 4: Fixed + Symmetric Problem (Maximum reduction)
# =============================================================================

print("=" * 70)
print("Example 4: FixedSymmetricProblem - Both constraints applied")
print("=" * 70)

form = FormDiagram.create_circular_radial_spaced(
    center=center,
    radius=radius,
    n_hoops=n_hoops,
    n_parallels=n_parallels,
    r_oculus=r_oculus,
)

for edge in form.edges_on_boundary():
    form.edge_attribute(edge, "_is_edge", False)

problem = FixedSymmetricProblem.from_formdiagram(form, method="QR", center=center, printout=False)

print(f"Vertices: {problem.n}")
print(f"Edges: {problem.m}")
print(f"Independent variables: {problem.k}")
print(f"Combined reduction: {100 * (1 - problem.k / problem.m):.1f}%")
print()

# =============================================================================
# Comparison Table
# =============================================================================

print("=" * 70)
print("Comparison Table")
print("=" * 70)
print(f"{'Problem Type':<25} {'Edges':<10} {'Variables':<12} {'Reduction':<12}")
print("-" * 70)

# Recreate all problems for fair comparison
form1 = FormDiagram.create_cross(x_span=x_span, y_span=y_span, n=n)
p1 = Problem.from_formdiagram(form1)
print(f"{'Problem':<25} {p1.m:<10} {p1.k:<12} {'-':<12}")

form2 = FormDiagram.create_cross(x_span=x_span, y_span=y_span, n=n)
p2 = FixedProblem.from_formdiagram(form2, method="QR")
print(f"{'FixedProblem':<25} {p2.m:<10} {p2.k:<12} {f'{100 * (1 - p2.k / p2.m):.1f}%':<12}")

form3 = FormDiagram.create_circular_radial_spaced(center=center, radius=radius, n_hoops=n_hoops, n_parallels=n_parallels, r_oculus=r_oculus)
for edge in form3.edges_on_boundary():
    form3.edge_attribute(edge, "_is_edge", False)
p3 = SymmetricProblem.from_formdiagram(form3, center=center)
print(f"{'SymmetricProblem':<25} {p3.m:<10} {p3.k:<12} {f'{100 * (1 - p3.k / p3.m):.1f}%':<12}")

form4 = FormDiagram.create_circular_radial_spaced(center=center, radius=radius, n_hoops=n_hoops, n_parallels=n_parallels, r_oculus=r_oculus)
for edge in form4.edges_on_boundary():
    form4.edge_attribute(edge, "_is_edge", False)
p4 = FixedSymmetricProblem.from_formdiagram(form4, method="QR", center=center)
print(f"{'FixedSymmetricProblem':<25} {p4.m:<10} {p4.k:<12} {f'{100 * (1 - p4.k / p4.m):.1f}%':<12}")

print()
print("=" * 70)
print("✓ All specialized problem classes demonstrated successfully!")
print("=" * 70)
