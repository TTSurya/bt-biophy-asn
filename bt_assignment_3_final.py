import numpy as np

# -------------------------------
# Data Structures
# -------------------------------
class Point:
    """Represents a 3D coordinate of a monomer."""
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z

# -------------------------------
# Input Data
# -------------------------------
sarw1 = [
    Point(6,1,0), Point(5,0,1), Point(6,1,2), Point(5,2,1), Point(4,3,0),
    Point(5,4,1), Point(6,5,2), Point(5,6,1), Point(4,5,0), Point(3,4,1), 
    Point(4,3,2), Point(5,2,3), Point(6,1,4), Point(5,2,5), Point(6,3,4), 
    Point(5,4,3), Point(4,3,4), Point(3,4,5), Point(2,5,4), Point(1,4,3), 
    Point(0,3,2), Point(1,2,3), Point(0,1,4), Point(1,0,3), Point(2,1,2)
]

sarw2 = [
    Point(4,1,2), Point(5,0,1), Point(6,1,2), Point(5,2,3), Point(4,3,4), 
    Point(3,2,5), Point(2,1,4), Point(3,0,3), Point(4,1,4), Point(3,2,3),
    Point(2,1,2), Point(3,2,1), Point(2,1,0), Point(1,2,1), Point(2,3,2), 
    Point(1,4,1), Point(2,5,0), Point(3,4,1), Point(4,3,0), Point(5,4,1), 
    Point(6,5,2), Point(5,6,3), Point(4,5,2), Point(5,6,1), Point(4,5,0)
]

# -------------------------------
# Parameters
# -------------------------------
complexity = 3   # distance cutoff for contactsp
n = len(sarw1)

# -------------------------------
# Helper Functions
# -------------------------------
def euclidean_distance(p1: Point, p2: Point) -> float:
    """Returns Euclidean distance between two monomers."""
    return np.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2)

def radius_of_gyration(chain: list[Point]) -> float:
    """Calculates the radius of gyration for a polymer chain."""
    coords = np.array([[p.x, p.y, p.z] for p in chain])
    center_of_mass = coords.mean(axis=0)
    rg_squared = np.mean(np.sum((coords - center_of_mass)**2, axis=1))
    return np.sqrt(rg_squared)

# -------------------------------
# Contact Calculations
# -------------------------------
is_correct = [False] * n
correct_contacts = 0

for i in range(n):
    dist = euclidean_distance(sarw1[i], sarw2[i])
    if dist**2 <= complexity:
        is_correct[i] = True
        correct_contacts += 1

incorrect_contacts = n - correct_contacts

# Trap contacts = correct contacts with neighbors also correct
trap_correct_contacts = sum(
    is_correct[i-1] and is_correct[i] and is_correct[i+1]
    for i in range(1, n-1)
)

# -------------------------------
# Radius of Gyration
# -------------------------------
rg1 = radius_of_gyration(sarw1)
rg2 = radius_of_gyration(sarw2)

# -------------------------------
# Output
# -------------------------------
print("--- Final Results --- 📊")
print(f"Correct Contacts     : {correct_contacts}")
print(f"Incorrect Contacts   : {incorrect_contacts}")
print(f"Trap Correct Contacts: {trap_correct_contacts}")
print("\n--- Radius of Gyration ---")
print(f"SARW 1 (Rg): {rg1:.4f}")
print(f"SARW 2 (Rg): {rg2:.4f}")
