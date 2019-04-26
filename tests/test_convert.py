from .context import readers
import os


class TestConvert:
    def __init__(self):
        pass

    @classmethod
    def setup_class(cls):
        testdir = os.path.abspath(os.path.dirname(__file__))
        test0 = os.path.join(testdir, "resources/test0.bngl")
        test1 = os.path.join(testdir, "resources/test1.bngl")
        cls.bm0 = readers.BNGLReader(test0).parse()
        cls.bm1 = readers.BNGLReader(test1).parse()

    @classmethod
    def teardown_class(cls):
        pass

    def test_bmodel0_read(self):
        assert len(self.bm0.parameters) == 8
        assert len(self.bm0.molecules) == 2
        assert len(self.bm0.initial_cond) == 2
        assert len(self.bm0.observables) == 3
        assert len(self.bm0.functions) == 1
        assert len(self.bm0.rules) == 1

    def test_bmodel1_read(self):
        assert len(self.bm1.parameters) == 29
        assert len(self.bm1.molecules) == 4
        assert len(self.bm1.initial_cond) == 4
        assert len(self.bm1.observables) == 7
        assert len(self.bm1.rules) == 19

    def test_bmodel0_convert(self):
        krules = [r for rule in self.bm0.rules for r in rule.convert()]
        assert len(krules) == 1

    def test_bmodel1_convert(self):
        krules = [r for rule in self.bm1.rules for r in rule.convert()]
        print('\n'.join([r.write_as_kappa() for r in krules]))
        assert len(krules) == 29
