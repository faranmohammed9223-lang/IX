import importlib
import unittest


class InstallSmokeTest(unittest.TestCase):
    def test_ixpy_import(self):
        module = importlib.import_module("ixpy")
        self.assertIsNotNone(module)

    def test_psdm_imports(self):
        psdm_module = importlib.import_module("PSDM.PSDM")
        functions_module = importlib.import_module("PSDM.PSDM_functions")
        self.assertTrue(hasattr(psdm_module, "PSDM"))
        self.assertTrue(hasattr(functions_module, "process_input_file"))


if __name__ == "__main__":
    unittest.main()
