"""Check pipe-flow results against known Moody-chart values."""
import unittest
from mechtools.pipe import colebrook, friction_factor, haaland, pressure_drop, reynolds


class FrictionFactor(unittest.TestCase):
    def test_laminar(self):
        self.assertAlmostEqual(friction_factor(1000), 0.064)

    def test_smooth_pipe_blasius_range(self):
        # Blasius: f = 0.316 Re^-0.25, good to ~2 % for smooth pipes at Re ≈ 1e4–1e5
        Re = 5e4
        self.assertAlmostEqual(friction_factor(Re, 0.0), 0.316 * Re**-0.25, delta=0.0006)

    def test_colebrook_reference_value(self):
        # Re = 1e5, ε/D = 0.001  ->  f ≈ 0.0222 (Moody chart / White, Fluid Mechanics)
        self.assertAlmostEqual(colebrook(1e5, 1e-3), 0.0222, delta=0.0002)

    def test_haaland_close_to_colebrook(self):
        for Re in (1e4, 1e5, 1e6):
            for rr in (0.0, 1e-4, 1e-2):
                c, h = colebrook(Re, rr), haaland(Re, rr)
                self.assertLess(abs(h - c) / c, 0.025)


class PressureDrop(unittest.TestCase):
    def test_water_in_steel_pipe(self):
        # 10 L/s of water through 100 m of 100 mm commercial steel pipe
        r = pressure_drop(0.010, 0.100, 100.0, roughness=4.5e-5)
        self.assertAlmostEqual(r.velocity, 1.273, places=3)
        self.assertEqual(r.regime, "turbulent")
        self.assertAlmostEqual(r.reynolds, reynolds(r.velocity, 0.1, 998.0, 1.002e-3))
        # hand calculation with f ≈ 0.0197 gives Δp ≈ 15.9 kPa
        self.assertAlmostEqual(r.pressure_drop / 1000, 15.9, delta=0.4)

    def test_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            friction_factor(0)


if __name__ == "__main__":
    unittest.main()
