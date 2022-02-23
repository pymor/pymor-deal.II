def test_demo_results(ndarrays_regression):
    from pymor_dealii.pymor.demo import run

    result, _, _, _ = run(plot_error=False)

    compare = ["errors", "basis_sizes", "rel_errors"]
    ndarrays_regression.check({k: v for k, v in result.items() if k in compare})
