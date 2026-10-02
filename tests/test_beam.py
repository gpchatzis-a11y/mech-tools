"""Check the numerical beam solver against closed-form textbook solutions."""
import unittest
from mechtools.beam import Beam, rectangular_I

E, I, L = 210e9, rectangular_I(0.05, 0.10), 4.0   # steel, 50 × 100 mm, 4 m


class SimplySupported(unittest.TestCase):
    def test_central_point_load(self):
        P = 10e3
        r = Beam(L, E, I).add_point_load(P, L / 2).solve()
        self.assertAlmostEqual(r.reactions["R_A"], P / 2)
        self.assertAlmostEqual(r.max_moment, P * L / 4, places=3)
        self.assertAlmostEqual(r.max_deflection, P * L**3 / (48 * E * I), delta=1e-7)
        self.assertAlmostEqual(r.max_deflection_at, L / 2, places=2)

    def test_uniform_load(self):
        w = 5e3
        r = Beam(L, E, I).add_uniform_load(w).solve()
        self.assertAlmostEqual(r.max_moment, w * L**2 / 8, places=2)
        self.assertAlmostEqual(r.max_deflection, 5 * w * L**4 / (384 * E * I), delta=1e-7)

    def test_offset_point_load_reactions(self):
        P, a = 12e3, 1.0
        r = Beam(L, E, I).add_point_load(P, a).solve()
        self.assertAlmostEqual(r.reactions["R_A"], P * (L - a) / L)
        self.assertAlmostEqual(r.reactions["R_B"], P * a / L)
        # max deflection of an offset load (Roark), b = distance to the nearer support:
        # P b (L²-b²)^1.5 / (9√3 L E I)
        b = min(a, L - a)
        expected = P * b * (L**2 - b**2) ** 1.5 / (9 * 3**0.5 * L * E * I)
        self.assertAlmostEqual(r.max_deflection, expected, delta=1e-7)


class Cantilever(unittest.TestCase):
    def test_end_point_load(self):
        P = 2e3
        r = Beam(L, E, I, support="cantilever").add_point_load(P, L).solve()
        self.assertAlmostEqual(r.reactions["M_A"], -P * L)
        self.assertAlmostEqual(r.max_deflection, P * L**3 / (3 * E * I), delta=1e-7)
        self.assertAlmostEqual(r.max_deflection_at, L)

    def test_uniform_load(self):
        w = 1e3
        r = Beam(L, E, I, support="cantilever").add_uniform_load(w).solve()
        self.assertAlmostEqual(r.max_deflection, w * L**4 / (8 * E * I), delta=1e-7)


class Validation(unittest.TestCase):
    def test_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            Beam(-1, E, I)
        with self.assertRaises(ValueError):
            Beam(L, E, I).add_point_load(1.0, L + 1)


if __name__ == "__main__":
    unittest.main()
