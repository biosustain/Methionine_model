.phony: figure_3 figure_4 sample_methionine sample_example_ode validate_methionine validate_example_ode build

MAIN_IDATA_6 = results/methionine/idata.nc
MAIN_IDATA_6_missing_ahcys = results/methionine_missing_ahcys/idata.nc
MAIN_IDATA_EXAMPLE_ODE = results/example_ode/idata.nc
MAIN_LAPLACE_EXAMPLE_ODE = results/example_ode/idata.nc
VALIDATE_6 = results/methionine_validation/user_input/config.toml
VALIDATE_6_missing_ahcys = results/methionine_missing_ahcys_validation/user_input/config.toml
FIGURE_3 = plots/figure_3.png
FIGURE_4 = plots/figure_4.png

$(MAIN_IDATA_6): | data/methionine/config.toml
	.venv/bin/maud sample data/methionine | tail -n 1 | xargs -I '{}' mv {}  results/methionine

$(MAIN_IDATA_6_missing_ahcys): | data/methionine_missing_ahcys/config.toml
	.venv/bin/maud sample data/methionine_missing_ahcys | tail -n 1 | xargs -I '{}' mv {}  results/methionine_missing_ahcys

$(MAIN_IDATA_EXAMPLE_ODE): | data/example_ode/config.toml
	.venv/bin/maud sample data/example_ode | tail -n 1 | xargs -I '{}' mv {} results/example_ode
	.venv/bin/maud laplace data/example_ode | tail -n 1 | xargs -I '{}' mv {} results/example_ode_laplace

$(VALIDATE_ODE): $(MAIN_IDATA_EXAMPLE_ODE)
	.venv/bin/maud predict results/example_ode | tail -n 1 | xargs -I "{}" mv {} results/example_ode_validation
	.venv/bin/maud predict results/example_ode_laplace | tail -n 1 | xargs -I "{}" mv {} results/example_ode_laplace_validation

$(VALIDATE_6): $(MAIN_IDATA_6)
	.venv/bin/maud predict results/methionine | tail -n 1 | xargs -I '{}' mv {}  results/methionine_validation

$(VALIDATE_6_missing_ahcys): $(MAIN_IDATA_6_missing_ahcys)
	.venv/bin/maud predict results/methionine_missing_ahcys | tail -n 1 | xargs -I '{}' mv {}  results/methionine_missing_ahcys_validation

$(FIGURE_2): $(MAIN_IDATA)

$(FIGURE_3): $(VALIDATE_ODE)

$(FIGURE_4): $(MAIN_IDATA_missing_ahcys)
	
figure_2: $(FIGURE_2)
	. .venv/bin/activate && ( python plot_figure_2.py)

figure_3: $(FIGURE_3)
	. .venv/bin/activate && ( python plot_figure_3.py)

figure_4: $(FIGURE_4)
	. .venv/bin/activate && ( python plot_figure_4.py)

sample_methionine: $(MAIN_IDATA_6) $(MAIN_IDATA_6_missing_ahcys)

sample_example_ode: $(MAIN_IDATA_EXAMPLE_ODE)

validate_methionine: $(VALIDATE_6) $(VALIDATE_6_missing_ahcys)

validate_example_ode: $(VALIDATE_ODE)

build: requirements.txt
	virtualenv .venv --prompt=maud
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install -r requirements.txt
