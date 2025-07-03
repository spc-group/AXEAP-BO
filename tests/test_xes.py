from axeap_bo import xes_calc
import pandas as pd
import numpy as np

# Test default spectral calculation for Ti1+
def test_ti_default():
    stick_x, stick_y = xes_calc.gen_sticks("Ti", "1+")
    m_data = pd.read_excel("./test_files/Ti1+_default_stick.xlsx")
    mstick_x = m_data["eV"].to_numpy()
    mstick_y = m_data["Intensity"].to_numpy()

    assert np.array_equal(stick_x, mstick_x)
    assert np.array_equal(stick_y, mstick_y)

# Test default spectral calculation for Cu3+

# Test custom multiplet calculations for Ti1+

# Test custom crystal field calculations for Ti1+

# Test custom energy range calculations for Ti1+

# Test 

# Test 