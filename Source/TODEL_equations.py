from sympy import *

# Forces
X, Y, U, V = symbols("X, Y, U, V")
POT = (X**2 + Y**2 + 2*X**2*Y - 2*Y**3/3)/2
KIN = (U**2 + V**2)/2
TOT = KIN + POT

DX = diff(POT, X)
DY = diff(POT, Y)

FX = -DX
FY = -DY

print(FX)
print(FY)

