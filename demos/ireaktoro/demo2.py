from reaktoro.interpreter.ireaktoro import interpret

input = """
ChemicalSystem:
    Database: supcrt98.xml
    AqueousPhase:
        Species: H2O(l) H+ OH- HCO3- CO2(aq) Ca+2 Cl- Na+ O2(aq) H2(aq)
    GaeousPhase:
        Species: H2O(g) CO2(g)
    MineralPhases: Calcite

Equilibrium State1:
    Temperature: 60 celsius
    Pressure: 150 bar
    Mixture:
        H2O: 1 kg
        NaCl: 1 umol
        HCl: 1 umol
        O2: 1 umol
        CaCO3: 10 mol

Equilibrium State2:
    Temperature: 60 celsius
    Pressure: 150 bar
    Mixture:
        H2O: 1 kg
        NaCl: 1 umol
        HCl: 1 mol
        O2: 1 umol
        CaCO3: 10 mol

EquilibriumPath:
    From: State1
    To: State2
    Plot 1:
        x: t
        y: m[Ca] pH
"""

interpret(input)
