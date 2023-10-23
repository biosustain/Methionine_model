.phony: figure_3 figure_4 sample validate build

MAIN_IDATA_6 = results/methionine/idata.nc
MAIN_IDATA_6_missing_ahcys = results/methionine_missing_ahcys/idata.nc
MAIN_IDATA_EXAMPLE_ODE = results/example_ode/idata.nc
MAIN_LAPLACE_EXAMPLE_ODE = results/example_ode/idata.nc
VALIDATE_6 = results/methionine_cycle_6_validation/user_input/config.toml
VALIDATE_6_missing_ahcys = results/methionine_cycle_6_missing_ahcys_validation/user_input/config.toml
FIGURE_3 = plots/figure_3.png
FIGURE_4 = plots/figure_4.png

$(MAIN_IDATA_6): data/methionine_cycle_6/config.toml
	.venv/bin/maud sample data/methionine_cycle | tail -n 1 | xargs -I '{}' mv {}  results/methionine_cycle_6

$(MAIN_IDATA_6_missing_ahcys): data/methionine_cycle_6_missing_ahcys/config.toml
	.venv/bin/maud sample data/methionine_cycle_missing_ahcys | tail -n 1 | xargs -I '{}' mv {}  results/methionine_cycle_6_missing_ahcys

$(MAIN_IDATA_EXAMPLE_ODE): data/example_ode/config.toml
	.venv/bin/maud sample data/example_ode | tail -n 1 | xargs -I '{}' mv {} results/example_ode
	.venv/bin/maud laplace data/example_ode | tail -n 1 | xargs -I '{}' mv {} results/example_ode_laplace

$(VALIDATE_ODE): $(MAIN_IDATA_EXAMPLE_ODE)
	.venv/bin/maud predict results/example_ode | tail -n 1 | xargs -I "{}" mv {} results/example_ode_validation
	.venv/bin/maud predict results/example_ode | tail -n 1 | xargs -I "{}" mv {} results/example_ode_laplace_validation

$(VALIDATE_6): $(MAIN_IDATA_6)
	.venv/bin/maud predict results/methionine_cycle_6 --oos_path="data/methionine_cycle_6" | tail -n 1 | xargs -I '{}' mv {}  results/methionine_cycle_6_validation

$(VALIDATE_6_missing_ahcys): $(MAIN_IDATA_6_missing_ahcys)
	.venv/bin/maud predict results/methionine_cycle_6_missing_ahcys --oos_path="data/methionine_cycle_6_missing_ahcys" | tail -n 1 | xargs -I '{}' mv {}  results/methionine_cycle_6_missing_ahcys_validation

$(FIGURE_3): $(MAIN_IDATA)
	.venv/bin/python plot_figure_3.py

$(FIGURE_4): $(VALIDATE_ODE)
	.venv/bin/python plot_figure_4.py
	
figure_3: $(FIGURE_3)

figure_4: $(FIGURE_4)

sample: $(MAIN_IDATA_6) $(MAIN_IDATA_6_missing_ahcys)

validate: $(VALIDATE_6) $(VALIDATE_6_missing_ahcys)

build:
	virtualenv .venv --prompt=maud
	.venv/bin/pip install maud-metabolic-models
	cp -rf .venv/lib/python*/site-packages/maud/data/example_inputs/methionine data
