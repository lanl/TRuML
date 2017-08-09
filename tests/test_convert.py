from .context import readers


class TestConvert:
    def __init__(self):
        pass

    @classmethod
    def setup_class(cls):
        cls.m0 = readers.BNGLReader("resources/test0.bngl").parse()
        cls.m1 = readers.BNGLReader("resources/test1.bngl").parse()

    @classmethod
    def teardown_class(cls):
        pass

    def test_model0(self):
        assert len(self.m0.parameters) == 8
        assert len(self.m0.molecules) == 2
        assert len(self.m0.initial_cond) == 2
        assert len(self.m0.observables) == 3
        assert len(self.m0.functions) == 1
        assert len(self.m0.rules) == 1

