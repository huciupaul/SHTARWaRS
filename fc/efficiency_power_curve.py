import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import pickle

# -------------------------------------------------------------------
# 1) "Digitized" (j, η) pairs for T=473 K (green‐square curve).
#
#    These 18 points were extracted by finding green‐square centroids
#    in the provided plot and mapping pixel‐coordinates back to data‐units.
#    You can tweak or add more points if you want a finer fit.
# -------------------------------------------------------------------
#
#   (j in A/m²,   η in %)

j_data = 20000 * np.array([0, 0.010121457177244596, 0.03238867068924164, 0.05566801447484513, 0.07692309771323576, 0.10829959951857593, 0.13765184077670323, 0.16700404342453545, 0.19028342582043403, 0.22064779735216775, 0.25607289747252376, 0.28036441014202873, 0.31072878167376244, 0.3562753389713631, 0.39271256936532545, 0.42206481062345275, 0.45242918215518646, 0.4898785428227554, 0.5283401109845208, 0.5647773413784832, 0.6022267020460521, 0.654858325700078, 0.692307686367647, 0.7297570470352158, 0.7601214185669496, 0.7894736598250769, 0.8198380313568105, 0.8572873920243795, 0.8876517635561131, 0.90283394932198, 0.9170040820348304, 0.940283464430729, 0.9554656501965959, 0.9696356284682661, 0.9777328250982981, 0.9848178142341331, 0.9929150108641647])
eta_data = 0.1 + 0.01 * np.array([77.5, 70.03027577703709, 63.85380474719604, 60.82612309628746, 58.888406377718425, 56.95069196908716, 55.37629843458981, 54.0441157362647, 53.196363950035185, 52.22750674571955, 51.016434085356124, 50.4108977551744, 49.563145968944895, 48.47318334654315, 47.625431560313636, 46.89878981204581, 46.29325348186409, 45.445501695634576, 44.59774990940506, 43.74999812317555, 42.902246336946035, 41.691173676582615, 40.84342651022866, 39.99567472399915, 39.14792293776963, 38.30017115154012, 37.452419365310604, 36.12024128686106, 34.78806320841151, 34.06141684026811, 33.33477509200028, 31.881486975589052, 30.42819885917783, 28.732699906594362, 27.40051720826925, 25.46280279963799, 22.67733660477722])

j_data_2 = 20000 * np.array([0, 0.05353534789173427, 0.09898989509774776, 0.12929292656842345, 0.17474747377443692, 0.24949496681502242, 0.3090909441202701, 0.35353540569017994, 0.4141414686315313, 0.4595960158375447, 0.5121212395608782, 0.5646464632842118, 0.6171717640721395, 0.6545454720601351, 0.6919192571127248, 0.7292929651007205, 0.7737374266706305, 0.7969697045593919, 0.8343434896119817, 0.8717171975999773, 0.9020202290706529])
P_data = np.array([0, 1000.0002825701673, 1685.1849026150162, 2129.629629629629, 2759.258835404007, 3740.7408820258233, 4481.481622766564, 5018.518235948349, 5685.185326470269, 6166.666525381584, 6703.703844988788, 7203.7032798484515, 7666.666242811414, 7962.962962962962, 8240.74074074074, 8499.999222932038, 8759.259117974176, 8851.851710566769, 8981.481481481482, 9037.036895751953, 8981.481481481482]) / 9037.036895751953


# Just to confirm we have 18 points:
# assert j_data.shape == eta_data.shape == (18,)
# -------------------------------------------------------------------
# 2) Build a cubic spline over [0, 20 000].
#
#    By default we will do a "natural" spline (second 
#    derivative = 0 at the two ends). If you want a different boundary
#    condition (e.g. "clamped"), you can pass bc_type="clamped" etc.
# -------------------------------------------------------------------
#
spline = CubicSpline(j_data, eta_data, bc_type='natural', extrapolate=True)
def eta_of_j(j):
    """
    Return the interpolated energy‐efficiency [%%] at current density j [A/m²]
    using a cubic‐spline based on the digitized T=473 K curve. 
    Returns NaN if j < 0 or j > 20 000 (no extrapolation).
    """
    j = np.atleast_1d(j)
    out = np.full_like(j, np.nan, dtype=float)
    mask = (j >= j_data.min()) & (j <= j_data.max())
    if np.any(mask):
        out[mask] = spline(j[mask])
    return out if out.shape != () else out.item()

spline2 = CubicSpline(j_data_2, P_data, bc_type='natural', extrapolate=True)
def P_of_j(j):
    """
    Return the interpolated power [W] at current density j [A/m²]
    using a cubic‐spline based on the digitized T=473 K curve. 
    Returns NaN if j < 0 or j > 20 000 (no extrapolation).
    """
    j = np.atleast_1d(j)
    out = np.full_like(j, np.nan, dtype=float)
    mask = (j >= j_data_2.min()) & (j <= j_data_2.max())
    if np.any(mask):
        out[mask] = spline2(j[mask])
    return out if out.shape != () else out.item()


# -------------------------------------------------------------------
# Combine splines to function eta_of_P
# --------------------------------------------------------------------
P_cropped = P_data[0:-1]
j_cropped = j_data_2[0:-1]
eta_cropped = eta_of_j(j_cropped)
print(eta_cropped)



spline_eta_of_P = CubicSpline(P_cropped, eta_cropped, bc_type='natural', extrapolate=True)
filepath = 'fc/spline_eta_of_P.pkl'
with open(filepath, 'wb') as f:
    pickle.dump(spline_eta_of_P, f, protocol=pickle.HIGHEST_PROTOCOL)


def eta_of_P(P):
    """
    Return the interpolated energy‐efficiency [%%] at power P [W]
    using a cubic‐spline based on the digitized T=473 K curve.
    Returns NaN if P < 0 or P > 9037.036895751953 (no extrapolation).
    """
    P = np.atleast_1d(P)
    out = np.full_like(P, np.nan, dtype=float)
    mask = (P >= P_cropped.min()) & (P <= P_cropped.max())
    if np.any(mask):
        out[mask] = spline_eta_of_P(P[mask])
    return out if out.shape != () else out.item()
    



# -------------------------------------------------------------------
# 3) (Optional) Quick plot to verify.
# -------------------------------------------------------------------
if __name__ == "__main__":
    # Plot original "dots" and the continuous spline
    P_dense = np.linspace(0, 1, 500)
    eta_dense = eta_of_P(P_dense)
    j_dense = np.linspace(0, 1, 500)
    eta_dense_j = eta_of_j(j_dense * 20000)
    plt.figure(figsize=(6,4))
    plt.plot(P_cropped, eta_cropped, 's', label="Digitized points (473 K)")
    plt.plot(P_dense, eta_dense, '-', label="Cubic‐spline, η(P)")
    plt.plot(j_dense , eta_dense_j, '-', label="Cubic‐spline, η(j)")
    plt.plot()
    plt.xlabel("Current density $j$ [A/m²]")
    plt.ylabel("Energy Efficiency η [%]")
    plt.title("Cubic‐Spline Interpolator for η at T=473 K")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()
    # # Example usage:
    # for test_j in [0, 1000, 5000, 10000, 15000, 20000]:
    #     print(f"η({test_j:.0f} A/m)  {eta_of_j(test_j):.2f} %")




