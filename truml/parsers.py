import pyparsing as pp


class KappaParser:

    name = pp.Word(pp.alphas, bodyChars=pp.alphanums + "_+-") | pp.Word("_", bodyChars=pp.alphanums + "_+-")

    lpar = pp.Suppress(pp.Literal("("))
    rpar = pp.Suppress(pp.Literal(")"))
    lbrace = pp.Suppress(pp.Literal("{"))
    rbrace = pp.Suppress(pp.Literal("}"))
    lbrak = pp.Suppress(pp.Literal("["))
    rbrak = pp.Suppress(pp.Literal("]"))
    comma = pp.Suppress(pp.Literal(","))
    uscore = pp.Literal("_")
    pound = pp.Literal("#")
    dot = pp.Literal(".")

    site_def = name.setResultsName('name') + \
           pp.Optional(lbrace + pp.Group(name + pp.ZeroOrMore(pp.Optional(comma) + name)).setResultsName('states') +
                       rbrace)
    agent_def = name.setResultsName('name') + lpar + \
            pp.Group(pp.Optional(site_def + pp.ZeroOrMore(pp.Optional(comma) + site_def))).setResultsName('sites') + rpar

    site_bond = lbrak + (uscore ^ pound ^ dot ^ pp.Word(pp.nums)) + rbrak
    site_state = lbrace + (name ^ dot) + rbrace
    site = name.setResultsName('name') + pp.Optional(site_state).setResultsName('state') + pp.Optional(site_bond).setResultsName('bond')

    agent = dot.setResultsName('name') ^ (name.setResultsName('name') + lpar + pp.Group(pp.Optional(site + pp.ZeroOrMore(pp.Optional(comma) + site))).setResultsName('sites') + rpar)

    def __init__(self):
        pass


class BNGLParser:

    name = pp.Word(pp.alphas, bodyChars=pp.alphanums + "_")

    lpar = pp.Suppress(pp.Literal("("))
    rpar = pp.Suppress(pp.Literal(")"))
    comma = pp.Suppress(pp.Literal(","))

    state_name = pp.Word(pp.alphanums, bodyChars=pp.alphanums + "_")
    state = pp.Suppress(pp.Literal("~")) + (state_name ^ pp.Literal("?"))
    bond_name = pp.Word(pp.nums + "+" + "?")
    bond = pp.Suppress(pp.Literal("!")) + bond_name

    site_def = name.setResultsName('name') + pp.ZeroOrMore(state).setResultsName('states')
    molecule_def = name.setResultsName('name') + lpar + pp.Optional(pp.delimitedList(site_def)).setResultsName('sites') + rpar

    site = name.setResultsName('name') + pp.Optional(state).setResultsName('state') + pp.Optional(bond).setResultsName('bond')
    molecule = name.setResultsName('name') + lpar + pp.Optional(pp.delimitedList(site)).setResultsName('sites') + rpar

    last_complex = molecule + pp.ZeroOrMore("." + molecule) + pp.White()
    last_complex.leaveWhitespace()

    def __init__(self):
        pass
