from .context import readers


class TestConvert:
    def __init__(self):
        pass

    @classmethod
    def setup_class(cls):
        cls.bm0 = readers.BNGLReader("tests/resources/test0.bngl").parse()
        cls.bm1 = readers.BNGLReader("tests/resources/test1.bngl").parse()

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
        mdefs = self.bm0.molecules
        krules = [r for rule in self.bm0.rules for r in rule.convert(mdefs, mdefs)]
        assert len(krules) == 1

    def test_bmodel1_convert(self):
        mdefs = self.bm1.molecules
        krules = [r for rule in self.bm1.rules for r in rule.convert(mdefs, mdefs)]
        print '\n'.join([x.write_as_kappa() for x in krules])
        assert len(krules) == 29
