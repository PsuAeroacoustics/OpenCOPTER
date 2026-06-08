# Loads real HART-II .tec measurement files and computes normalized L2 error.
# This will  tests the file-reading + L2 validation logic using actual HART-II data files.
#not FINAL-prot2
import sys
import numpy as np


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


def load_tec_file(file_path):
    numeric_rows = []

    with open(file_path, "r") as file:
        for line in file:
            parts = line.split()

            try:
                row = [float(value) for value in parts]
                numeric_rows.append(row)
            except ValueError:
                # to skip header lines like TITLE, VARIABLES, ZONE
                pass

    return np.array(numeric_rows)


def main():
    baseline_path = "oc_fly/example/hart_ii/mn-contour-meas.tec"
    new_path = "oc_fly/example/hart_ii/mv-contour-meas.tec"

    baseline_data = load_tec_file(baseline_path)
    new_data = load_tec_file(new_path)

    print("Baseline file:", baseline_path)
    print("New file:", new_path)
    print("Baseline shape:", baseline_data.shape)
    print("New shape:", new_data.shape)

    baseline_signal = baseline_data[:, 5]
    new_signal = new_data[:, 5]

    validator = L2Validator(threshold=0.15)
    result = validator.validate(new_signal, baseline_signal)

    print("Compared column: OASPL1to244")
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

