import math
import unittest
from pathlib import Path


AMREX_ROOT = Path(__file__).resolve().parents[2]
HEADER = AMREX_ROOT / "Src" / "AmrCore" / "AMReX_Interpolater.H"
SOURCE = AMREX_ROOT / "Src" / "AmrCore" / "AMReX_Interpolater.cpp"
LINEAR_WEIGHTS = (3.0 / 8.0, 3.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0)
EPS = 1.0e-12


def symmetric_z_tau(beta):
    return 0.5 * abs((beta[0] + beta[1]) - (beta[2] + beta[3]))


def js_weights(beta):
    alpha = [
        linear / (indicator + EPS) ** 2
        for linear, indicator in zip(LINEAR_WEIGHTS, beta)
    ]
    total = sum(alpha)
    return tuple(value / total for value in alpha)


class TestHWENOWeights(unittest.TestCase):
    def test_symmetric_z_tau_is_reflection_invariant(self):
        beta = (0.2, 1.7, 0.4, 2.3)
        reflected = (beta[1], beta[0], beta[3], beta[2])
        self.assertEqual(symmetric_z_tau(beta), symmetric_z_tau(reflected))

    def test_js_weights_use_inverse_squared_smoothness(self):
        beta = (0.2, 0.5, 1.0, 2.0)
        weights = js_weights(beta)
        self.assertTrue(all(math.isfinite(weight) for weight in weights))
        self.assertAlmostEqual(sum(weights), 1.0, places=15)
        self.assertGreater(weights[0] / LINEAR_WEIGHTS[0], weights[1] / LINEAR_WEIGHTS[1])
        self.assertGreater(weights[1] / LINEAR_WEIGHTS[1], weights[2] / LINEAR_WEIGHTS[2])

    def test_header_keeps_symmetric_z_and_exposes_js_weights(self):
        header = HEADER.read_text(encoding="utf-8")
        self.assertIn(
            "(beta[0] + beta[1]) - (beta[2] + beta[3])", header
        )
        self.assertIn("const Real tau = Real(0.5) *", header)
        self.assertIn("static GpuArray<Real,4> JSWeights", header)
        self.assertIn("d[k] * inv_beta * inv_beta", header)

    def test_prolong_uses_js_while_restriction_keeps_z(self):
        header = HEADER.read_text(encoding="utf-8")
        weighted_start = header.index("static GpuArray<Real,4> WeightedCubic")
        weighted_end = header.index("static Real hweno_child_coordinate", weighted_start)
        weighted_cubic = header[weighted_start:weighted_end]
        self.assertIn("const auto omega = JSWeights(beta);", weighted_cubic)
        self.assertNotIn("const auto omega = ZWeights(beta);", weighted_cubic)

        source = SOURCE.read_text(encoding="utf-8")
        restrict_start = source.index("HermiteWENO2D::hweno_restrict_y")
        restrict_end = source.index(
            "HermiteWENO2D::configure_amr_transfer_pp", restrict_start
        )
        restriction = source[restrict_start:restrict_end]
        self.assertEqual(
            restriction.count("const auto omega = ZWeights(beta);"), 2
        )
        self.assertNotIn("JSWeights(beta)", restriction)


if __name__ == "__main__":
    unittest.main()
