SHELL = /bin/bash

install:
	@./check_moog.sh
	@mkdir -p linelist
	@mkdir -p spectra
	@tar -zxvf models/apogee_kurucz.tar.gz
	@mv apogee_kurucz models
	@tar -zxvf models/marcs.tar.gz
	@mv marcs models
	@echo "Atmosphere models installed in dir: models"
	@echo "Installing dependencies..."
	@pip install -r requirements.txt
	@conda install -c anaconda wxpython=3.0.0.0
	@echo "Dependencies installed"
	@echo ""
	@echo "MOOGme is successfully installed!"
	@echo "Type"
	@echo "    python FASMA.py"
	@echo "to start. Happy spectroscopying :)"
