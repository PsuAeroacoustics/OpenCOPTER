import sys
import argparse
import numpy as np
from scipy.io import loadmat


class ValidationResult:
    def __init__(self, error_value, threshold, passed):
        self.error_value = error_value
        self.threshold = threshold
        self.passed = passed


class L2Validator:
    def __init__(self, threshold=0.15):
        self.threshold = threshold

    def l2_norm(self, signal):
        return float(np.sqrt(np.sum(signal * signal)))

    def normalized_l2_error(self, new_signal, baseline_signal):

        new_signal = np.ravel(new_signal)
        baseline_signal = np.ravel(baseline_signal)

        n = min(len(new_signal), len(baseline_signal))

        new_signal = new_signal[:n]
        baseline_signal = baseline_signal[:n]

        difference = new_signal - baseline_signal

        numerator = self.l2_norm(difference)
        denominator = self.l2_norm(baseline_signal) + 1e-12

        return numerator / denominator

    def validate(self, new_signal, baseline_signal):
        error = self.normalized_l2_error(new_signal, baseline_signal)
        passed = error <= self.threshold
        return ValidationResult(error, self.threshold, passed)


def load_mat_signal(mat_path, key):

    data = loadmat(mat_path)

    if key not in data:
        raise KeyError(f"{key} not found in {mat_path}")

    return np.asarray(data[key], dtype=np.float64)


def main():

    parser = argparse.ArgumentParser(
        description="Validate OpenCOPTER MAT outputs using normalized L2 error"
    )

    parser.add_argument("--baseline", required=True)
    parser.add_argument("--new", required=True)
    parser.add_argument("--key", required=True)
    parser.add_argument("--threshold", type=float, default=0.15)

    args = parser.parse_args()

    baseline_signal = load_mat_signal(args.baseline, args.key)
    new_signal = load_mat_signal(args.new, args.key)

    print("Baseline file:", args.baseline)
    print("New file:", args.new)
    print("Key:", args.key)

    print("Baseline shape:", baseline_signal.shape)
    print("New shape:", new_signal.shape)

    validator = L2Validator(threshold=args.threshold)

    result = validator.validate(new_signal, baseline_signal)

    print("Normalized L2 error:", result.error_value)
    print("Threshold:", result.threshold)

    if result.passed:
        print("PASS")
        return 0
    else:
        print("FAIL")
        return 1


if __name__ == "__main__":
    sys.exit(main())
