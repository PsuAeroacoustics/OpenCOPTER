# validation script for HART-II
# this will compare a new predicted signal to a baseline signal using normalized L2 error
#not the FINAL VERSION-prot1

import sys
import numpy as np




class ValidationResult:
    # object to hold the final validation result
    def __init__(self, error_value, threshold, passed):
        self.error_value = error_value
        self.threshold = threshold
        self.passed = passed


class L2Validator:
    # object that computes L2 error and decides pass/fail
    def __init__(self, threshold=0.02):
        self.threshold = threshold

    def l2_norm(self, signal):
        # L2 norm = sqrt(sum(signal_i^2))
        squared_values = signal * signal
        total = np.sum(squared_values)
        return float(np.sqrt(total))

    def normalized_l2_error(self, new_signal, baseline_signal):
        # making  both signals the same length before comparing
        n = min(len(new_signal), len(baseline_signal))

        new_signal = new_signal[:n]
        baseline_signal = baseline_signal[:n]

        difference = new_signal - baseline_signal

        numerator = self.l2_norm(difference)
        denominator = self.l2_norm(baseline_signal) + 1e-12

        return numerator / denominator

    def validate(self, new_signal, baseline_signal):
        error = self.normalized_l2_error(new_signal, baseline_signal)

        if error <= self.threshold:
            passed = True
        else:
            passed = False

        return ValidationResult(error, self.threshold, passed)


class HartIIDemo:
    # demo class using placeholder arrays until real HART-II files can be loaded
    def __init__(self):
        self.baseline_signal = np.array([1.1, 2.0, 2.9, 4.2])
        self.new_signal = np.array([1.0, 2.0, 3.0, 4.0])

        self.validator = L2Validator(threshold=0.15)

    def run(self):
        result = self.validator.validate(self.new_signal, self.baseline_signal)

        print("HART-II validation demo")
        print("Normalized L2 error:", result.error_value)
        print("Threshold:", result.threshold)

        if result.passed:
            print("PASS")
            return 0
        else:
            print("FAIL")
            return 1


def main():
    demo = HartIIDemo()
    return demo.run()


if __name__ == "__main__":
    sys.exit(main())

