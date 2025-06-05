import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
# -------------------------------------------------------------------
# 1) "Digitized" (j, η) pairs for T=473 K (green‐square curve).
#
#    These 18 points were extracted by finding green‐square centroids
#    in the provided plot and mapping pixel‐coordinates back to data‐units.
#    You can tweak or add more points if you want a finer fit.
# -------------------------------------------------------------------
#
#   (j in A/m²,   η in %)

j_data = 20000 * np.array([0.010121457177244596, 0.03238867068924164, 0.05566801447484513, 0.07692309771323576, 0.10829959951857593, 0.13765184077670323, 0.16700404342453545, 0.19028342582043403, 0.22064779735216775, 0.25607289747252376, 0.28036441014202873, 0.31072878167376244, 0.3562753389713631, 0.39271256936532545, 0.42206481062345275, 0.45242918215518646, 0.4898785428227554, 0.5283401109845208, 0.5647773413784832, 0.6022267020460521, 0.654858325700078, 0.692307686367647, 0.7297570470352158, 0.7601214185669496, 0.7894736598250769, 0.8198380313568105, 0.8572873920243795, 0.8876517635561131, 0.90283394932198, 0.9170040820348304, 0.940283464430729, 0.9554656501965959, 0.9696356284682661, 0.9777328250982981, 0.9848178142341331, 0.9929150108641647])
eta_data = ([70.03027577703709, 63.85380474719604, 60.82612309628746, 58.888406377718425, 56.95069196908716, 55.37629843458981, 54.0441157362647, 53.196363950035185, 52.22750674571955, 51.016434085356124, 50.4108977551744, 49.563145968944895, 48.47318334654315, 47.625431560313636, 46.89878981204581, 46.29325348186409, 45.445501695634576, 44.59774990940506, 43.74999812317555, 42.902246336946035, 41.691173676582615, 40.84342651022866, 39.99567472399915, 39.14792293776963, 38.30017115154012, 37.452419365310604, 36.12024128686106, 34.78806320841151, 34.06141684026811, 33.33477509200028, 31.881486975589052, 30.42819885917783, 28.732699906594362, 27.40051720826925, 25.46280279963799, 22.67733660477722])
# Just to confirm we have 18 points:
assert j_data.shape == eta_data.shape == (18,)
# -------------------------------------------------------------------
# 2) Build a cubic spline over [0, 20 000].
#
#    By default we will do a "natural" spline (second 
#    derivative = 0 at the two ends). If you want a different boundary
#    condition (e.g. "clamped"), you can pass bc_type="clamped" etc.
# -------------------------------------------------------------------
#
spline = CubicSpline(j_data, eta_data, bc_type='natural', extrapolate=False)
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
# -------------------------------------------------------------------
# 3) (Optional) Quick plot to verify.
# -------------------------------------------------------------------
if __name__ == "__main__":
    # Plot original "dots" and the continuous spline
    j_dense = np.linspace(0, 20000, 500)
    eta_dense = eta_of_j(j_dense)
    plt.figure(figsize=(6,4))
    plt.plot(j_data, eta_data, 's', label="Digitized points (473 K)")
    plt.plot(j_dense, eta_dense, '-', label="Cubic‐spline, η(j)")
    plt.xlabel("Current density $j$ [A/m²]")
    plt.ylabel("Energy Efficiency ηₜ꜀ [%]")
    plt.title("Cubic‐Spline Interpolator for ηₜ꜀ at T=473 K")
    plt.xlim(0, 20000)
    plt.ylim(0, 75)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()
    # Example usage:
    for test_j in [0, 1000, 5000, 10000, 15000, 20000]:
        print(f"η({test_j:.0f} A/m²) ≃ {eta_of_j(test_j):.2f} %")